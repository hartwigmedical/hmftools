package com.hartwig.hmftools.svtools.rna_expression;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;

public class GeneReadData
{
    public final EnsemblGeneData GeneData;

    private final List<RegionReadData> mRegionReadData;

    // summary results
    private int mTotalReadCount;
    private List<TranscriptResults> mTranscriptResults;

    public GeneReadData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;

        mRegionReadData = Lists.newArrayList();

        mTranscriptResults = Lists.newArrayList();
        mTotalReadCount = 0;
    }

    public final List<RegionReadData> getRegionReadData() { return mRegionReadData; }

    public void addRegionReadData(final RegionReadData regionReadData)
    {
        if(hasRegionData(regionReadData.Region.start(), regionReadData.Region.end()))
            return;

        mRegionReadData.add(regionReadData);
    }

    public boolean hasRegionData(long posStart, long posEnd)
    {
        return mRegionReadData.stream().anyMatch(x -> x.Region.start() == posStart && x.Region.end() == posEnd);
    }

    public RegionReadData findRegionData(long posStart, long posEnd)
    {
        return mRegionReadData.stream()
                .filter(x -> x.Region.start() == posStart && x.Region.end() == posEnd)
                .findFirst().orElse(null);
    }

    public int totalReadCount() { return mTotalReadCount; }
    public void addReadCount() { ++mTotalReadCount; }
    public void setTotalReadCount(int count) { mTotalReadCount = count; }

    public List<TranscriptResults> getTranscriptResults() { return mTranscriptResults; }
}
