package com.hartwig.hmftools.sage.dedup;

import static com.hartwig.hmftools.sage.dedup.VariantDeduper.longerContainsShorter;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.DEDUP_MNV_FILTER;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.DEDUP_SNV_MNV_FILTER;

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

            index = nextIndex;
        }
    }

    private static void dedupMnvPair(final SageVariant first, final SageVariant second)
    {
        int firstRcLength = first.readContext().indexedBases().coreLength();
        int secondRcLength = second.readContext().indexedBases().coreLength();

        if(firstRcLength > secondRcLength)
        {
            second.filters().add(DEDUP_MNV_FILTER);
        }
        else if(firstRcLength < secondRcLength)
        {
            first.filters().add(DEDUP_MNV_FILTER);
        }
        else
        {
            // equal so use tumor quality
            if(first.totalQuality() > second.totalQuality())
                second.filters().add(DEDUP_MNV_FILTER);
            else
                first.filters().add(DEDUP_MNV_FILTER);
        }
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
            snv.filters().add(DEDUP_SNV_MNV_FILTER);
            return true;
        }

        return false;
    }
}
