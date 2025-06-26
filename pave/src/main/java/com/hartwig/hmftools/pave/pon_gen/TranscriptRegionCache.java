package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.Collections;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class TranscriptRegionCache
{
    private final List<BaseRegion> mExonicRegions;
    private int mCurrentIndex;

    private static final int SPLICE_BUFFER = 10;

    public TranscriptRegionCache(final List<BaseRegion> regions)
    {
        mExonicRegions = regions;
        mCurrentIndex = 0;
    }

    public static TranscriptRegionCache from(final EnsemblDataCache ensemblDataCache, final ChrBaseRegion region)
    {
        List<BaseRegion> exonicRegions = Lists.newArrayList();

        if(ensemblDataCache == null)
            return new TranscriptRegionCache(Collections.emptyList());

        List<GeneData> geneDataList = ensemblDataCache.getChrGeneDataMap().get(region.Chromosome);

        if(geneDataList != null && !geneDataList.isEmpty())
        {
            for(GeneData geneData : geneDataList)
            {
                if(!region.overlaps(geneData.Chromosome, geneData.GeneStart, geneData.GeneEnd))
                    continue;

                TranscriptData transcriptData = ensemblDataCache.getCanonicalTranscriptData(geneData.GeneId);

                if(transcriptData != null)
                {
                    for(ExonData exonData : transcriptData.exons())
                    {
                        // add a splice region buffer
                        if(positionsOverlap(region.start(), region.end(), exonData.Start - SPLICE_BUFFER, exonData.End + SPLICE_BUFFER))
                        {
                            exonicRegions.add(new BaseRegion(exonData.Start - SPLICE_BUFFER, exonData.End + SPLICE_BUFFER));
                        }
                    }
                }
            }
        }

        return new TranscriptRegionCache(exonicRegions);
    }

    public void resetSearch() { mCurrentIndex = 0; }
    public int entryCount() { return mExonicRegions.size(); }

    public boolean inOrNearExonicRegion(int position)
    {
        if(mExonicRegions.isEmpty())
            return false;

        // cannot search at earlier positions than the current index
        if(mCurrentIndex < mExonicRegions.size() && position < mExonicRegions.get(mCurrentIndex).start())
            return false;

        boolean matched = false;

        for(; mCurrentIndex < mExonicRegions.size(); ++mCurrentIndex)
        {
            BaseRegion exonRegion = mExonicRegions.get(mCurrentIndex);

            if(exonRegion.end() < position)
                continue;

            if(exonRegion.start() > position)
            {
                // move index back
                if(mCurrentIndex > 0)
                    --mCurrentIndex;

                return false;
            }

            matched = true;
            break;
        }

        return matched;
    }

    @VisibleForTesting
    public int currentIndex() { return mCurrentIndex; }
}
