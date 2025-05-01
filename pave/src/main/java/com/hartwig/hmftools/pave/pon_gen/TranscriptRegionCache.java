package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.List;

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

    public TranscriptRegionCache(final EnsemblDataCache ensemblDataCache, final ChrBaseRegion region)
    {
        mExonicRegions = Lists.newArrayList();
        mCurrentIndex = 0;

        if(ensemblDataCache == null)
            return;

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
                            mExonicRegions.add(new BaseRegion(exonData.Start - SPLICE_BUFFER, exonData.End + SPLICE_BUFFER));
                        }
                    }
                }
            }
        }
    }

    public void resetSearch() { mCurrentIndex = 0; }
    public int entryCount() { return mExonicRegions.size(); }

    public boolean inOrNearExonicRegion(int position)
    {
        if(mExonicRegions == null || mExonicRegions.isEmpty())
            return false;

        int firstPosMatchIndex = -1;
        boolean matched = false;

        for(; mCurrentIndex < mExonicRegions.size(); ++mCurrentIndex)
        {
            BaseRegion exonRegion = mExonicRegions.get(mCurrentIndex);

            if(exonRegion.end() < position)
                continue;

            if(exonRegion.start() > position)
                break;

            if(firstPosMatchIndex == -1)
                firstPosMatchIndex = mCurrentIndex;

            matched = true;
            break;
        }

        // move the index back to the prior position or the first at this position
        if(firstPosMatchIndex >= 0)
            mCurrentIndex = firstPosMatchIndex;
        else if(mCurrentIndex > 0)
            --mCurrentIndex;

        return matched;
    }
}
