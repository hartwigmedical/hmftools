package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Nucleotides.isValidDnaBase;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_INVALID_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_HP_LEN;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractTpValues;
import static com.hartwig.hmftools.sage.seqtech.Homopolymer.getHomopolymers;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.Microhomology;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.samtools.SAMRecord;

public final class UltimaUtils
{
    protected static final int MAX_HOMOPOLYMER = 15;

    protected static final byte INVALID_BASE = -1;
    private static final byte TP_ZERO_BASE_QUAL = 0;

    // equivalent to adding their logs, ie P(combined) = 10^(-qual/10) + 10^(-qual/10), combined qual = -10 * log10(P(combined))
    protected static final int HALF_PHRED_SCORE_SCALING = 3;

    protected static final UltimaQualRecalibration BQR_CACHE = new UltimaQualRecalibration();

    public static void loadBqrCache(final String filename)
    {
        BQR_CACHE.loadRecalibrationFile(filename);
    }

    public static void setMaxRawQual(final byte qual) { BQR_CACHE.setMaxRawQual(qual); }
    public static byte maxRawQual() { return BQR_CACHE.maxRawQual(); }

    public static byte safeQualLookup(final byte[] qualityArray, int index)
    {
        if(index < 0 || index >= qualityArray.length)
            return ULTIMA_INVALID_QUAL;

        return qualityArray[index];
    }

    public static byte calcTpBaseQual(final SAMRecord record, int indexStart, int indexEnd, int tpSearchValue)
    {
        if(indexStart < 0 || indexStart >= record.getReadBases().length)
            return ULTIMA_INVALID_QUAL;

        byte[] tpValues = extractTpValues(record);

        byte qualValue1 = -1;
        byte qualValue2 = -1;

        for(int i = indexStart; i <= min(indexEnd, tpValues.length - 1); ++i)
        {
            if(tpValues[i] == tpSearchValue)
            {
                if(qualValue1 < 0)
                {
                    qualValue1 = record.getBaseQualities()[i];
                }
                else
                {
                    qualValue2 = record.getBaseQualities()[i];
                    break;
                }
            }
        }

        int homopolymerLength = indexEnd - indexStart + 1;
        char homopolymerBase = (char)record.getReadBases()[indexStart];

        // check quals vs BQR
        qualValue1 = BQR_CACHE.calcTpRecalibratedQual(qualValue1, homopolymerLength, homopolymerBase, false);
        qualValue2 = BQR_CACHE.calcTpRecalibratedQual(qualValue2, homopolymerLength, homopolymerBase, false);

        if(qualValue1 < 0)
        {
            // if bases aren't found, use middle base(s) as: min(40, Q + 6 * abs(required_tp – T))
            int middleIndex = indexStart + (indexEnd - indexStart) / 2;
            final byte[] baseQualities = record.getBaseQualities();
            if(middleIndex < 0 || middleIndex >= baseQualities.length)
            {
                return ULTIMA_INVALID_QUAL;
            }

            qualValue1 = baseQualities[middleIndex];
            byte tpValue = tpValues[middleIndex];

            if(tpValue == TP_ZERO_BASE_QUAL)
                return BQR_CACHE.calcTpRecalibratedQual(qualValue1, homopolymerLength, homopolymerBase, true);;

            qualValue1 = BQR_CACHE.calcTpRecalibratedQual(qualValue1, homopolymerLength, homopolymerBase, false);
            byte bqrValue = BQR_CACHE.calcTpRecalibratedQual(BQR_CACHE.maxRawQual(), homopolymerLength, homopolymerBase, false);

            return (byte)min(bqrValue, qualValue1 + HALF_PHRED_SCORE_SCALING * 2 * abs(tpSearchValue - tpValue));
        }

        if(qualValue2 < 0)
            return qualValue1;

        return (byte)(qualValue1 - HALF_PHRED_SCORE_SCALING);
    }

    protected static int findHomopolymerLength(final byte[] refBases, final byte compareBase, int startIndex, boolean searchUp)
    {
        int repeatCount = 0;

        int i = startIndex;

        while(i >= 0 && i < refBases.length)
        {
            if(refBases[i] != compareBase)
                return repeatCount;

            ++repeatCount;

            if(searchUp)
                ++i;
            else
                --i;
        }

        return repeatCount;
    }

    public static List<Integer> coreHomopolymerLengths(final VariantReadContext readContext)
    {
        List<Homopolymer> homopolymers = getHomopolymers(readContext.ReadBases, readContext.CoreIndexStart, readContext.CoreIndexEnd);
        return homopolymers.stream().map(x -> x.Length).collect(Collectors.toList());
    }

    @VisibleForTesting
    public static boolean isCleanSnv(final VariantReadContext readContext)
    {
        if(!readContext.variant().isSNV())
        {
            return false;
        }

        if(readContext.coreLength() != readContext.refBasesBytes().length)
        {
            return false;
        }

        for(int i = readContext.CoreIndexStart; i <= readContext.CoreIndexEnd; ++i)
        {
            if(i == readContext.VarIndex)
            {
                continue;  // the SNV base
            }

            byte readBase = readContext.readBasesBytes()[i];
            byte refBase = readContext.refBasesBytes()[i - readContext.CoreIndexStart];
            if(readBase != refBase)
            {
                return false;
            }
        }

        return true;
    }

    @VisibleForTesting
    public static boolean isMsiIndelOfType(final SimpleVariant variant, final List<String> units)
    {
        if(!variant.isIndel())
            return false;

        String indelBases = variant.isInsert() ? variant.alt().substring(1) : variant.ref().substring(1);
        for(String unit : units)
        {
            int numRepeats = indelBases.length() / unit.length();
            String expectedIndelBases = unit.repeat(numRepeats);
            if(indelBases.equals(expectedIndelBases))
                return true;
        }

        return false;
    }

    @VisibleForTesting
    public static boolean isAdjacentToLongHomopolymer(final SAMRecord read, int varIndexInRead, int longLength)
    {
        byte[] readBases = read.getReadBases();
        boolean roomOnRight = varIndexInRead + longLength < read.getReadBases().length;
        if(roomOnRight)
        {
            int startIndex = varIndexInRead + 1;
            int endIndex = varIndexInRead + longLength;
            if(isDnaBaseHomopolymer(readBases, startIndex, endIndex))
                return true;
        }

        boolean roomOnLeft = varIndexInRead >= longLength;
        if(roomOnLeft)
        {
            int startIndex = varIndexInRead - longLength;
            int endIndex = varIndexInRead - 1;
            if(isDnaBaseHomopolymer(readBases, startIndex, endIndex))
                return true;
        }

        return false;
    }

    private static boolean isDnaBaseHomopolymer(final byte[] bases, int startIndex, int endIndex)
    {
        if(!isValidDnaBase(bases[startIndex]))
            return false;

        for(int i = startIndex + 1; i <= endIndex; i++)
        {
            if(bases[i - 1] != bases[i])
                return false;
        }

        return true;
    }

    private static final List<String> SINGLE_HOMOPOLYMERS = List.of("A", "C", "G", "T");
    private static final List<String> DI_NUC_HOMOPOLYMERS = List.of("TA", "AT");

    public static boolean ultimaLongRepeatFilter(
            final SimpleVariant variant, final SAMRecord read, int varIndexInRead, final Microhomology homology)
    {
        if(isMsiIndelOfType(variant, SINGLE_HOMOPOLYMERS) && homology != null && homology.Length >= ULTIMA_MAX_HP_LEN - 5)
            return true; // Ultima cannot call variants of this type

        if(isMsiIndelOfType(variant, DI_NUC_HOMOPOLYMERS) && homology != null && homology.Length >= ULTIMA_MAX_HP_LEN)
            return true; // Ultima cannot call variants of this type

        if(isAdjacentToLongHomopolymer(read, varIndexInRead, ULTIMA_MAX_HP_LEN))
            return true;  // Ultima gets confused near variants of this type

        return false;
    }

    // variant filters
    public static boolean belowExpectedHpQuals(final ReadContextCounter primaryTumor)
    {
        if(!primaryTumor.isIndel())
            return false;

        UltimaVariantData ultimaData = primaryTumor.ultimaData();
        final List<Integer> homopolymerLengths = ultimaData.homopolymerLengths();

        // TODO: make constants and generally improve
        if(primaryTumor.isLongIndel() && Collections.max(homopolymerLengths) < 5)
            return false;

        for(int i = 0; i < homopolymerLengths.size(); i++)
        {
            int length = homopolymerLengths.get(i);

            double avgQual = ultimaData.homopolymerAvgQuals().get(i);

            if(length == 1 && avgQual < 24)
                return true;
            else if(length == 2 && avgQual < 22)
                return true;
            else if(length == 3 && avgQual < 18)
                return true;
            else if(length == 4 && avgQual < 18)
                return true;
            else if(length == 5 && avgQual < 16)
                return true;
            else if(length == 6 && avgQual < 14)
                return true;
            else if(length == 7 && avgQual < 12)
                return true;
            else if(avgQual < 10)
                return true;
            else if(length >= 15)
                return true;
        }
        return false;
    }

    public static boolean belowExpectedT0Quals(final ReadContextCounter primaryTumor, final boolean nearbyVariant)
    {
        // TODO: make constants
        int threshold = nearbyVariant ? 28 : 18;
        return Collections.min(primaryTumor.ultimaData().t0AvgQuals()) < threshold;
    }
}
