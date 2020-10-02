package com.hartwig.hmftools.isofox.common;

import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionsOverlap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionReadData.regionExists;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;
import static com.hartwig.hmftools.isofox.common.RnaUtils.deriveCommonRegions;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GeneReadData
{
    public final EnsemblGeneData GeneData;

    private final List<RegionReadData> mExonRegions; // set of unique exons ie with differing start and end positions

    private final List<TranscriptData> mTranscripts;
    private boolean mHasUnsplicedRegions;

    public GeneReadData(final EnsemblGeneData geneData)
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

        // check for any region which doesn't overlap any others
        for(RegionReadData region : mExonRegions)
        {
            if(mExonRegions.stream()
                    .filter(x -> x != region)
                    .noneMatch(x -> positionsOverlap(region.start(), region.end(), x.start(), x.end())))
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

    public static List<GeneReadData> createGeneReadData(final List<EnsemblGeneData> geneDataList, final EnsemblDataCache geneTransCache)
    {
        final List<GeneReadData> geneReadDataList = Lists.newArrayList();

        for(EnsemblGeneData geneData : geneDataList)
        {
            GeneReadData geneReadData = new GeneReadData(geneData);

            List<TranscriptData> transDataList = Lists.newArrayList(geneTransCache.getTranscripts(geneData.GeneId));

            if(transDataList.isEmpty())
            {
                ISF_LOGGER.warn("no transcripts found for gene({}:{})", geneData.GeneId, geneData.GeneName);
                continue;
            }

            geneReadData.setTranscripts(transDataList);
            geneReadDataList.add(geneReadData);
        }

        return geneReadDataList;
    }

}
