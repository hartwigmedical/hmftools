package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.codon.Nucleotides.isValidDnaBase;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.HALF_PHRED_SCORE_SCALING;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_INVALID_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_HP_LEN;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractTpValues;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.filter.FilterConfig.ULTIMA_CANDIDATE_HIGH_BQ_REPEAT_MIN;
import static com.hartwig.hmftools.sage.seqtech.UltimaQualModelBuilder.canSkipRealignedModels;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import static com.hartwig.hmftools.sage.filter.VariantFilters.STRAND_BIAS_CALCS;
import com.hartwig.hmftools.sage.common.Microhomology;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.samtools.SAMRecord;

public final class UltimaUtils
{
    protected static final int MAX_HOMOPOLYMER = 15;

    protected static final byte INVALID_BASE = -1;
    private static final byte TP_ZERO_BASE_QUAL = 0;

    private static final byte T0_EXPECTED_QUAL_THRESHOLD_SIMPLE = 18;
    private static final byte T0_EXPECTED_QUAL_THRESHOLD_COMPLEX = 24;

    protected static final UltimaQualRecalibration BQR_CACHE = new UltimaQualRecalibration();

    public static void loadBqrCache(final String filename)
    {
        List<String> lines = null;

        if(filename != null)
        {
            try
            {
                lines = Files.readAllLines(new File(filename).toPath());
            }
            catch(Exception e)
            {
                SG_LOGGER.error("failed to read Ultima base-qual recalibration file({}): {}", filename, e.toString());
            }
        }
        else
        {
            InputStream inputStream = UltimaUtils.class.getResourceAsStream("/seqtech/ultima_qual_recalibration.tsv");
            lines = new BufferedReader(new InputStreamReader(inputStream)).lines().collect(toList());
        }

        BQR_CACHE.loadRecalibrationData(lines);
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
        boolean isReverse = record.getReadNegativeStrandFlag();

        // check quals vs BQR
        qualValue1 = BQR_CACHE.calcTpRecalibratedQual(qualValue1, homopolymerLength, homopolymerBase, false, isReverse);
        qualValue2 = BQR_CACHE.calcTpRecalibratedQual(qualValue2, homopolymerLength, homopolymerBase, false, isReverse);

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
                return BQR_CACHE.calcTpRecalibratedQual(qualValue1, homopolymerLength, homopolymerBase, true, isReverse);

            qualValue1 = BQR_CACHE.calcTpRecalibratedQual(qualValue1, homopolymerLength, homopolymerBase, false, isReverse);

            byte bqrValue = BQR_CACHE.calcTpRecalibratedQual(
                    BQR_CACHE.maxRawQual(), homopolymerLength, homopolymerBase, false, isReverse);

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

    public static boolean hasLongHomopolymerInContext(final SimpleVariant variant, final SAMRecord read, final int varReadIndex)
    {
        int standardLength = MIN_CORE_DISTANCE + DEFAULT_FLANK_LENGTH;
        int coreStartIndex = max(varReadIndex - standardLength + (variant.isIndel() ? 1 : 0), 0);
        int coreEndIndex = min(varReadIndex + (variant.altLength() - 1) + standardLength, read.getReadBases().length - 1);

        byte lastBase = read.getReadBases()[coreStartIndex];
        int repeatCount = 1;
        for(int i = coreStartIndex + 1; i <= coreEndIndex; ++i)
        {
            byte nextBase = read.getReadBases()[i];
            if(nextBase == lastBase)
            {
                ++repeatCount;

                if(repeatCount >= ULTIMA_CANDIDATE_HIGH_BQ_REPEAT_MIN)
                    return true;
            }
            else
            {
                lastBase = nextBase;
                repeatCount = 1;
            }
        }

        return false;
    }

    public static boolean lowQualInReadCore(final SimpleVariant variant, int variantReadIndex, final List<Integer> lowQualIndices)
    {
        if(lowQualIndices == null || lowQualIndices.isEmpty())
            return false;

        int coreIndexStart = variantReadIndex - MIN_CORE_DISTANCE + (variant.isIndel() ? 1 : 0);
        int coreIndexEnd = variantReadIndex + variant.alt().length() - 1 + MIN_CORE_DISTANCE;

        for(int index = coreIndexStart; index <= coreIndexEnd; ++index)
        {
            if(lowQualIndices.contains(index))
                return true;
        }

        return false;
    }

    public static boolean lowQualInReadContextCore(int variantReadIndex, final List<Integer> lowQualIndices, final VariantReadContext readContext)
    {
        if(lowQualIndices == null || lowQualIndices.isEmpty())
            return false;

        // exclude from checking if the variant has a long homopolymer
        if(readContext.MaxRepeat != null && readContext.MaxRepeat.repeatLength() == 1
        && readContext.MaxRepeat.Count >= ULTIMA_CANDIDATE_HIGH_BQ_REPEAT_MIN)
        {
            return false;
        }

        // check for a sufficiently long single-base repeat anywhere in the read bases
        byte lastBase = readContext.ReadBases[0];
        int repeatCount = 1;
        for(int i = 1; i < readContext.ReadBases.length; ++i)
        {
            byte nextBase = readContext.ReadBases[i];
            if(nextBase == lastBase)
            {
                ++repeatCount;

                if(repeatCount >= ULTIMA_CANDIDATE_HIGH_BQ_REPEAT_MIN)
                    return false;
            }
            else
            {
                lastBase = nextBase;
                repeatCount = 1;
            }
        }

        int coreIndexStart = variantReadIndex - readContext.leftCoreLength();
        int coreIndexEnd = variantReadIndex + readContext.rightCoreLength();

        for(int index = coreIndexStart; index <= coreIndexEnd; ++index)
        {
            if(lowQualIndices.contains(index))
                return true;
        }

        return false;
    }

    private static final List<String> SINGLE_HOMOPOLYMERS = List.of("A", "C", "G", "T");
    private static final List<String> DI_NUCLEOTIDES = List.of("TA", "AT");

    public static boolean ultimaLongRepeatFilter(
            final SimpleVariant variant, final SAMRecord read, int varIndexInRead, final Microhomology homology)
    {
        if(isMsiIndelOfType(variant, SINGLE_HOMOPOLYMERS) && homology != null && homology.Length >= ULTIMA_MAX_HP_LEN - 5)
            return true; // Ultima cannot call variants of this type

        if(isMsiIndelOfType(variant, DI_NUCLEOTIDES) && homology != null && homology.Length >= ULTIMA_MAX_HP_LEN)
            return true; // Ultima cannot call variants of this type

        if(isAdjacentToLongHomopolymer(read, varIndexInRead, ULTIMA_MAX_HP_LEN))
            return true;  // Ultima gets confused near variants of this type

        return false;
    }

    // variant filters
    private static final Map<Integer,Double> MIN_HP_QUAL_MAP = Maps.newHashMap();

    static
    {
        MIN_HP_QUAL_MAP.put(1, 21.5);
        MIN_HP_QUAL_MAP.put(2, 20.0);
        MIN_HP_QUAL_MAP.put(3, 18.5);
        MIN_HP_QUAL_MAP.put(4, 17.0);
        MIN_HP_QUAL_MAP.put(5, 15.5);
        MIN_HP_QUAL_MAP.put(6, 14.0);
        MIN_HP_QUAL_MAP.put(7, 12.0);
    }

    public static boolean belowExpectedHpQuals(final ReadContextCounter primaryTumor)
    {
        if(primaryTumor.isSnv() && !primaryTumor.readContext().hasIndelInCore())
            return false;

        UltimaVariantData ultimaData = primaryTumor.ultimaData();
        List<Integer> homopolymerLengths = ultimaData.paddedHomopolymerLengths();

        // TODO: make constants and generally improve
        if(primaryTumor.isLongIndel() && Collections.max(homopolymerLengths) < 5)
            return false;

        for(int i = 0; i < homopolymerLengths.size(); i++)
        {
            int length = homopolymerLengths.get(i);

            if(length >= 15)
                return true;

            double threshold = MIN_HP_QUAL_MAP.getOrDefault(length, 10.0);
            if(!primaryTumor.isIndel())
                threshold -= 2.5;

            double avgQual = ultimaData.homopolymerAvgQuals().get(i);

            if(avgQual < threshold)
                return true;
        }

        return false;
    }

    public static boolean belowExpectedT0Quals(final ReadContextCounter primaryTumor)
    {
        boolean isSimple = canSkipRealignedModels(primaryTumor.readContext());

        if(primaryTumor.variant().isMNV())
            isSimple = !STRAND_BIAS_CALCS.allOneSide(
                    primaryTumor.readStrandBiasAlt()) && !primaryTumor.ultimaData().getQualModels().mOriginalQualModel.canCompute();

        int threshold = isSimple ? T0_EXPECTED_QUAL_THRESHOLD_SIMPLE : T0_EXPECTED_QUAL_THRESHOLD_COMPLEX;

        return Collections.min(primaryTumor.ultimaData().t0AvgQuals()) < threshold;
    }
}
