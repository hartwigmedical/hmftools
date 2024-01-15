package com.hartwig.hmftools.sage.dedup;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.DEDUP_INDEL_FILTER_OLD;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.SageVariant;

public final class DedupIndelOld
{
    /*
    - make a direct comparison of variants in the same LPS where the CORE at least partially overlap
      and where at least 1 of the variants must be an INDEL.
      For each case where this occurs:
        - Calculate the read context and flank of each variant excluding the variant
        - If the CORE of one variant is fully explained by the CORE+FLANKS of the other, then filter as DEDUP whichever has a shorter CORE
        or if identical, the lowest quality.
    */

    public static void dedupIndelsOld(final List<SageVariant> variants)
    {
        List<SageVariant> candidates = variants.stream()
                .filter(x -> x.isPassing())
                .filter(x -> x.hasLocalPhaseSets())
                .collect(Collectors.toList());

        for(int i = 0; i < candidates.size() - 1; ++i)
        {
            SageVariant firstVariant = candidates.get(i);
            int firstCoreEnd = firstVariant.readContext().indexedBases().corePositionEnd();

            for(int j = i + 1; j < candidates.size(); ++j)
            {
                SageVariant nextVariant = candidates.get(j);

                // cores must overlap or a delete contains another variant
                int nextCoreStart = nextVariant.readContext().indexedBases().corePositionStart();
                boolean firstDeletesSecond = firstVariant.isDelete() && deleteContainsVariant(firstVariant, nextVariant);

                if(nextCoreStart > firstCoreEnd && !firstDeletesSecond)
                    break;

                if(!validPair(firstVariant, nextVariant))
                    continue;

                dedupPair(firstVariant, nextVariant);
            }
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

    private static void dedupPair(final SageVariant first, final SageVariant second)
    {
        // calculate the read context and flank of each variant excluding the variant
        // if the CORE of one variant is fully explained by the CORE+FLANKS of the other,
        // or if an INDEL coexist at the same base as another variant and shares an LPS, or a DEL overlaps any other variant which shares an LPS
        // then the variant with the shorter read context should be filtered as DEDUP

        // filter as DEDUP whichever has a shorter CORE or if identical, the lowest quality

        final CoreFlanksRefBases firstRefBases = extractCoreFlanksLessAlt(first);
        final CoreFlanksRefBases secondRefBases = extractCoreFlanksLessAlt(second);

        boolean applyFilter = false;

        if(secondRefBases.CorePlusFlanks.contains(firstRefBases.Core)
        || firstRefBases.CorePlusFlanks.contains(secondRefBases.Core))
        {
            applyFilter = true;
        }
        else if(first.position() == second.position())
        {
            applyFilter = true;
        }
        else if(first.isDelete() && deleteContainsVariant(first, second))
        {
            applyFilter = true;
        }

        if(!applyFilter)
            return;

        // first check for matching indels, and then take the left-aligned one
        if(first.alt().length() == second.alt().length() && first.ref().length() == second.ref().length())
        {
            second.filters().add(DEDUP_INDEL_FILTER_OLD);
            return;
        }

        int firstRcLength = first.readContext().indexedBases().length();
        int secondRcLength = second.readContext().indexedBases().length();

        if(secondRcLength < firstRcLength)
        {
            second.filters().add(DEDUP_INDEL_FILTER_OLD);
        }
        else if(secondRcLength > firstRcLength)
        {
            first.filters().add(DEDUP_INDEL_FILTER_OLD);
        }
        else
        {
            if(first.totalQuality() >= second.totalQuality())
                second.filters().add(DEDUP_INDEL_FILTER_OLD);
            else
                first.filters().add(DEDUP_INDEL_FILTER_OLD);
        }
    }

    private static boolean deleteContainsVariant(final SageVariant del, final SageVariant variant)
    {
        int delStart = del.position() + 1;
        int delEnd = del.position() + del.ref().length() - 1;

        return variant.position() >= delStart && variant.position() + variant.ref().length() - 1 <= delEnd;
    }

    private static CoreFlanksRefBases extractCoreFlanksLessAlt(final SageVariant variant)
    {
        final IndexedBases indexedBases = variant.readContext().indexedBases();

        try
        {
            String coreString = indexedBases.coreString();
            int altIndex = indexedBases.Index - indexedBases.LeftCoreIndex;
            int postAltIndex = altIndex + variant.alt().length();
            String coreRefString = coreString.substring(0, altIndex) + variant.ref();

            if(postAltIndex <= coreString.length())
                coreRefString += coreString.substring(postAltIndex);

            String coreFlankRefBases = indexedBases.leftFlankString() + coreRefString + indexedBases.rightFlankString();
            return new CoreFlanksRefBases(coreRefString, coreFlankRefBases);
        }
        catch(Exception e)
        {
            SG_LOGGER.error("var({}) error extracting core and flanks: {}",
                    variant, indexedBases);
            return new CoreFlanksRefBases("", "");
        }
    }

    static final class CoreFlanksRefBases
    {
        public final String Core;
        public final String CorePlusFlanks;

        public CoreFlanksRefBases(final String core, final String corePlusFlanks)
        {
            Core = core;
            CorePlusFlanks = corePlusFlanks;
        }
    }
}
