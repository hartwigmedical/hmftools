package com.hartwig.hmftools.isofox.common;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionReadData.regionExists;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.RnaUtils.deriveCommonRegions;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;

public class GeneReadData
{
    public final com.hartwig.hmftools.common.gene.GeneData GeneData;

    private final List<RegionReadData> mExonRegions; // set of unique exons ie with differing start and end positions

    private final List<TranscriptData> mTranscripts;
    private boolean mHasUnsplicedRegions;

    public GeneReadData(final com.hartwig.hmftools.common.gene.GeneData geneData)
    {
        GeneData = geneData;

        mExonRegions = Lists.newArrayList();
        mTranscripts = Lists.newArrayList();
        mHasUnsplicedRegions = false;
    }

    public String name() { return GeneData.GeneName;}

    public final List<TranscriptData> getTranscripts() { return mTranscripts; }
    public void setTranscripts(final List<TranscriptData> transDataList)
    {
        mTranscripts.addAll(transDataList);
    }

    public final List<RegionReadData> getExonRegions() { return mExonRegions; }

    public void addExonRegion(final RegionReadData region)
    {
        if(!regionExists(mExonRegions, region.start(), region.end()))
            mExonRegions.add(region);
    }

    public boolean hasUnsplicedRegions() { return mHasUnsplicedRegions; }

    public void setHasUnsplicedRegions()
    {
        mHasUnsplicedRegions = false;

        if(mExonRegions.isEmpty())
            return;

        int transPosMin = mExonRegions.stream().mapToInt(x -> x.start()).min().orElse(0);
        int transPosMax = mExonRegions.stream().mapToInt(x -> x.end()).max().orElse(0);

        for(RegionReadData region : mExonRegions)
        {
            // if any base directly before or after this exon isn't covered by another region, then it is unspliced
            if(region.start() > transPosMin
            && mExonRegions.stream().noneMatch(x -> positionWithin(region.start() - 1, x.start(), x.end())))
            {
                mHasUnsplicedRegions = true;
                return;
            }

            if(region.end() < transPosMax
            && mExonRegions.stream().noneMatch(x -> positionWithin(region.end() + 1, x.start(), x.end())))
            {
                mHasUnsplicedRegions = true;
                return;
            }
        }
    }

    public static void generateCommonExonicRegions(final List<RegionReadData> regions, final List<int[]> allCommonRegions)
    {
        if(regions.isEmpty())
            return;

        List<int[]> commonRegions = Lists.newArrayList(new int[] {regions.get(0).start(), regions.get(0).end()});

        for(int i = 1; i < regions.size(); ++i)
        {
            List<int[]> nextRegion = Lists.newArrayList(new int[] {regions.get(i).start(), regions.get(i).end()});
            commonRegions = deriveCommonRegions(commonRegions, nextRegion);
        }

        allCommonRegions.addAll(commonRegions);
    }

    public int calcExonicRegionLength()
    {
        final List<int[]> commonExonicRegions = Lists.newArrayList();
        generateCommonExonicRegions(mExonRegions, commonExonicRegions);
        return commonExonicRegions.stream().mapToInt(x -> x[SE_END] - x[SE_START]).sum();
    }

    public List<RegionReadData> findOverlappingRegions(final ReadRecord read)
    {
        return ReadRecord.findOverlappingRegions(mExonRegions, read);
    }

    public String toString()
    {
        return String.format("%s:%s location(%s:%d -> %d) trans(%d)",
                GeneData.GeneId, GeneData.GeneName, GeneData.Chromosome, GeneData.GeneStart, GeneData.GeneEnd, mTranscripts.size());
    }

    public static List<GeneReadData> createGeneReadData(final List<com.hartwig.hmftools.common.gene.GeneData> geneDataList, final EnsemblDataCache geneTransCache)
    {
        final List<GeneReadData> geneReadDataList = Lists.newArrayList();

        for(com.hartwig.hmftools.common.gene.GeneData geneData : geneDataList)
        {
            GeneReadData geneReadData = new GeneReadData(geneData);

            final List<TranscriptData> geneTranscripts = geneTransCache.getTranscripts(geneData.GeneId);

            if(geneTranscripts == null || geneTranscripts.isEmpty())
            {
                ISF_LOGGER.warn("no transcripts found for gene({}:{} chr={})", geneData.GeneId, geneData.GeneName, geneData.Chromosome);
                continue;
            }

            geneReadData.setTranscripts(geneTranscripts);
            geneReadDataList.add(geneReadData);
        }

        return geneReadDataList;
    }



}
