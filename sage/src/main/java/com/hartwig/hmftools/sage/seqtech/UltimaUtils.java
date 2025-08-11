package com.hartwig.hmftools.sage.seqtech;

import static com.hartwig.hmftools.common.codon.Nucleotides.isValidDnaBase;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.CYCLE_BASES;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_BOOSTED_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_INVALID_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL_T0;
import static com.hartwig.hmftools.sage.seqtech.Homopolymer.getHomopolymers;

import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.Microhomology;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public final class UltimaUtils
{
    protected static final int MAX_HOMOPOLYMER = 15;

    protected static final byte INVALID_BASE = -1;

    private static final int LONG_LENGTH = 20;

    public static byte safeQualLookup(final byte[] qualityArray, int index)
    {
        if(index < 0 || index >= qualityArray.length)
            return ULTIMA_INVALID_QUAL;

        byte value = qualityArray[index];
        if(value == ULTIMA_BOOSTED_QUAL)
            return ULTIMA_MAX_QUAL_T0;

        return value;
    }

    protected static int findHomopolymerLength(final byte[] refBases, final byte compareBase, int startIndex, boolean searchUp)
    {
        // byte repeatBase = refBases[startIndex];
        int repeatCount = 0;

        // int i = startIndex + (searchUp ? 1 : -1);
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

    public static boolean ultimaLongRepeatFilter(
            final SimpleVariant variant, final SAMRecord read, int varIndexInRead, final Microhomology homology)
    {
        if(isMsiIndelOfType(variant, List.of("A", "C", "G", "T")) && homology != null && homology.Length >= LONG_LENGTH - 5)
            return true; // Ultima cannot call variants of this type

        if(isMsiIndelOfType(variant, List.of("TA", "AT")) && homology != null && homology.Length >= LONG_LENGTH)
            return true; // Ultima cannot call variants of this type

        if(isAdjacentToLongHomopolymer(read, varIndexInRead, LONG_LENGTH))
            return true;  // Ultima gets confused near variants of this type

        return false;
    }
}
