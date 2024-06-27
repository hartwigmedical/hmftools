package com.hartwig.hmftools.sage.dedup;

import static com.hartwig.hmftools.sage.dedup.VariantDeduper.longerContainsShorter;
import static com.hartwig.hmftools.sage.filter.SoftFilter.DEDUP_MNV;
import static com.hartwig.hmftools.sage.filter.SoftFilter.DEDUP_SNV_MNV;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.sage.common.SageVariant;

public final class DedupSnvMnv
{
    public static void dedupMnvOverlaps(final List<SageVariant> variants)
    {
        List<SageVariant> mnvs = variants.stream()
                .filter(x -> x.isMnv())
                .filter(x -> x.isPassing())
                .filter(x -> x.hasLocalPhaseSets())
                .collect(Collectors.toList());;

        int index = 0;
        while(index < mnvs.size() - 1)
        {
            SageVariant firstVariant = mnvs.get(index);

            // look for an overlapping MNV
            int maxPos = firstVariant.position() + firstVariant.alt().length() - 1;
            int nextIndex = index + 1;

            while(nextIndex < mnvs.size())
            {
                SageVariant nextVariant = mnvs.get(nextIndex);

                if(nextVariant.position() > maxPos)
                    break;

                if(!nextVariant.hasMatchingLps(firstVariant.localPhaseSets()))
                {
                    ++nextIndex;
                    continue;
                }

                // filter the overlaps
                DedupSnvMnv.dedupMnvPair(firstVariant, nextVariant);

                if(!firstVariant.isPassing())
                    break;

                ++nextIndex;
            }

            ++index;
        }
    }

    private static void dedupMnvPair(final SageVariant first, final SageVariant second)
    {
        // First any MNV which has 2 changed bases and overlaps with an MNV with 3 changed bases is filtered.
        // If both have 2 changed bases then the least compact MNV is filtered
        // If both have the same number of bases and changed bases then the lowest qual MNV is filtered
        int firstChangedBases = mnvChangedBaseCount(first);
        int secondChangedBases = mnvChangedBaseCount(second);

        if(firstChangedBases > secondChangedBases)
        {
            second.filters().add(DEDUP_MNV);
        }
        else if(secondChangedBases > firstChangedBases)
        {
            first.filters().add(DEDUP_MNV);
        }
        else
        {
            if(first.ref().length() < second.ref().length())
            {
                second.filters().add(DEDUP_MNV);
            }
            else if(second.ref().length() < first.ref().length())
            {
                first.filters().add(DEDUP_MNV);
            }
            else
            {
                // equal so use tumor quality
                if(first.totalQuality() > second.totalQuality())
                    second.filters().add(DEDUP_MNV);
                else
                    first.filters().add(DEDUP_MNV);
            }
        }
    }

    private static int mnvChangedBaseCount(final SageVariant mnv)
    {
        int changed = 0;
        for(int i = 0; i < mnv.ref().length(); ++i)
        {
            if(mnv.ref().charAt(i) != mnv.alt().charAt(i))
                ++changed;
        }

        return changed;
    }

    public static void dedupMnvSnvs(final List<SageVariant> variants)
    {
        List<SageVariant> candidates = variants.stream()
                .filter(x -> !x.isIndel())
                .filter(x -> x.isPassing())
                .filter(x -> x.hasLocalPhaseSets())
                .collect(Collectors.toList());;

        int index = 0;
        while(index < candidates.size() - 1)
        {
            SageVariant firstVariant = candidates.get(index);

            // look for an overlapping MNV
            int maxPos = firstVariant.position() + firstVariant.alt().length() - 1;
            int nextIndex = index + 1;

            while(nextIndex < candidates.size())
            {
                SageVariant nextVariant = candidates.get(nextIndex);

                if(nextVariant.position() > maxPos) // must overlap
                    break;

                if(!nextVariant.hasMatchingLps(firstVariant.localPhaseSets()))
                {
                    ++nextIndex;
                    continue;
                }

                // filter the overlaps
                DedupSnvMnv.dedupMnvSnvPair(firstVariant, nextVariant);

                if(!firstVariant.isPassing())
                    break;

                ++nextIndex;
            }

            ++index;
        }
    }

    private static boolean dedupMnvSnvPair(final SageVariant first, final SageVariant second)
    {
        if(first.isSnv() == second.isSnv())
            return false;

        SageVariant snv = first.isSnv() ? first : second;
        SageVariant mnv = first.isMnv() ? first : second;

        if(longerContainsShorter(snv, mnv))
        {
            snv.filters().add(DEDUP_SNV_MNV);
            return true;
        }

        return false;
    }
}
