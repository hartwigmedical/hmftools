package com.hartwig.hmftools.sage.dedup;

import static com.hartwig.hmftools.sage.vcf.VariantVCF.DEDUP_INDEL_FILTER;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.DEDUP_MNV_FILTER;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.DEDUP_SNV_MNV_FILTER;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.SageVariant;

public final class DedupIndel
{
    /*
    - replaces test for local realign set, DedupIndel and DedupRealign
    - Replace the concept of ‘LRS’ with direct comparison of variants in the same LPS where the CORE at least partially overlap
      and where at least 1 of the variants must be an INDEL.
      For each case where this occurs:
        - Calculate the read context and flank of each variant excluding the variant
        - If the CORE of one variant is fully explained by the CORE+FLANKS of the other, then filter as DEDUP whichever has a shorter CORE
        or if identical, the lowest quality.
    */

    public static void dedupIndels(final List<SageVariant> variants)
    {
        List<SageVariant> candidates = variants.stream()
                .filter(x -> x.isPassing())
                .filter(x -> x.hasLocalPhaseSets())
                .collect(Collectors.toList());;

        int index = 0;
        while(index < candidates.size() - 1)
        {
            SageVariant firstVariant = candidates.get(index);

            // look for an overlapping MNV
            int firstCoreEnd = firstVariant.readContext().indexedBases().corePositionEnd();

            int nextIndex = index + 1;

            while(nextIndex < candidates.size())
            {
                SageVariant nextVariant = candidates.get(nextIndex);

                // cores must overlap
                int nextCoreStart = nextVariant.readContext().indexedBases().corePositionStart();
                if(nextCoreStart > firstCoreEnd)
                    break;

                if(!validPair(firstVariant, nextVariant))
                {
                    ++nextIndex;
                    continue;
                }

                dedupPair(firstVariant, nextVariant);
                ++nextIndex;
            }

            index = nextIndex;
        }

    }

    private static boolean validPair(final SageVariant first, final SageVariant second)
    {
        if(!first.isPassing() || !second.isPassing())
            return false;

        if(!first.hasMatchingLps(second.localPhaseSets()))
            return false;

        if(!first.isIndel() && !second.isIndel())
            return false;

        return true;
    }

    private static boolean dedupPair(final SageVariant first, final SageVariant second)
    {
        // calculate the read context and flank of each variant excluding the variant
        // if the CORE of one variant is fully explained by the CORE+FLANKS of the other, then
        // filter as DEDUP whichever has a shorter CORE or if identical, the lowest quality

        final String[] firstRefBases = extractCoreFlanksLessAlt(first);
        final String[] secondRefBases = extractCoreFlanksLessAlt(second);

        if(!secondRefBases[CORE_FLANKS_STR].contains(firstRefBases[CORE_STR])
        && !firstRefBases[CORE_FLANKS_STR].contains(secondRefBases[CORE_STR]))
        {
            return false;
        }

        int firstRcLength = first.readContext().indexedBases().length();
        int secondRcLength = second.readContext().indexedBases().length();

        if(secondRcLength < firstRcLength)
        {
            second.filters().add(DEDUP_INDEL_FILTER);
        }
        else if(secondRcLength > firstRcLength)
        {
            first.filters().add(DEDUP_INDEL_FILTER);
        }
        else
        {
            if(first.totalQuality() > second.totalQuality())
                second.filters().add(DEDUP_MNV_FILTER);
            else
                first.filters().add(DEDUP_MNV_FILTER);
        }

        return false;
    }

    private static final int CORE_STR = 0;
    private static final int CORE_FLANKS_STR = 1;

    private static String[] extractCoreFlanksLessAlt(final SageVariant variant)
    {
        final IndexedBases indexedBases = variant.readContext().indexedBases();
        String coreString = indexedBases.coreString();
        int altIndex = indexedBases.Index - indexedBases.LeftCoreIndex;
        int postAltIndex = altIndex + variant.alt().length();
        String coreRefString = coreString.substring(0, altIndex) + variant.ref() + coreString.substring(postAltIndex);
        String coreFlankRefBases = indexedBases.leftFlankString() + coreRefString + indexedBases.rightFlankString();
        return new String[] { coreRefString, coreFlankRefBases };
    }

}
