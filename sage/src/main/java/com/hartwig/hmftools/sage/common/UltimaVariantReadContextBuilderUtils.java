package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.SimpleVariant;

import htsjdk.samtools.SAMRecord;

public class UltimaVariantReadContextBuilderUtils
{
    private static final int LONG_LENGTH = 20;

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
        if(baseIndex(bases[startIndex]) < 0)
            return false;

        for(int i = startIndex + 1; i <= endIndex; i++)
        {
            if(bases[i - 1] != bases[i])
                return false;
        }

        return true;
    }

    public static boolean ultimaLongRepeatFilter(final SimpleVariant variant, final SAMRecord read, int varIndexInRead,
            final Microhomology homology)
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
