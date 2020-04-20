package com.hartwig.hmftools.isofox.common;

import static com.hartwig.hmftools.isofox.common.RegionReadData.regionExists;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;
import static com.hartwig.hmftools.isofox.common.RnaUtils.deriveCommonRegions;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GeneReadData
{
    public final EnsemblGeneData GeneData;

    private final List<RegionReadData> mExonRegions; // set of unique exons ie with differing start and end positions

    private final List<TranscriptData> mTranscripts;

    public GeneReadData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;

        mExonRegions = Lists.newArrayList();
        mTranscripts = Lists.newArrayList();
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
        if(!regionExists(mExonRegions, region.PosStart, region.PosEnd))
            mExonRegions.add(region);
    }

    public static void generateCommonExonicRegions(final List<RegionReadData> regions, final List<long[]> allCommonRegions)
    {
        if(regions.isEmpty())
            return;

        List<long[]> commonRegions = Lists.newArrayList(new long[] {regions.get(0).start(), regions.get(0).end()});

        for(int i = 1; i < regions.size(); ++i)
        {
            List<long[]> nextRegion = Lists.newArrayList(new long[] {regions.get(i).start(), regions.get(i).end()});
            commonRegions = deriveCommonRegions(commonRegions, nextRegion);
        }

        allCommonRegions.addAll(commonRegions);
    }

    public long calcExonicRegionLength()
    {
        final List<long[]> commonExonicRegions = Lists.newArrayList();
        generateCommonExonicRegions(mExonRegions, commonExonicRegions);
        return commonExonicRegions.stream().mapToLong(x -> x[SE_END] - x[SE_START]).sum();
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

}
