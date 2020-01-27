package com.hartwig.hmftools.svtools.rna_expression;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import htsjdk.samtools.SAMRecord;

public class RegionReadData
{
    public final GenomeRegion Region;

    public final String Id;

    private String mRefBases;
    private int[] mRefBasesMatched;
    private int mMatchingReads;
    private Map<RegionReadData,Integer> mLinkedRegions; // where a read covers this region and another next to it

    public RegionReadData(final GenomeRegion region, final String id)
    {
        Region = region;
        Id = id;

        mRefBasesMatched = new int[(int)region.bases()+1];
        mRefBases = "";
        mLinkedRegions = Maps.newHashMap();
        mMatchingReads = 0;
    }

    public String chromosome() { return Region.chromosome(); }
    public long start() { return Region.start(); }
    public long end() { return Region.end(); }

    public int[] refBasesMatched() { return mRefBasesMatched; }
    public void addMatchedRead() { ++mMatchingReads; }
    public int matchedReadCount() { return mMatchingReads; }

    public final Map<RegionReadData, Integer> getLinkedRegions() { return mLinkedRegions; }

    public void addLinkedRegion(final RegionReadData region)
    {
        Integer linkCount = mLinkedRegions.get(region);

        if(linkCount != null)
            mLinkedRegions.put(region, linkCount + 1);
        else
            mLinkedRegions.put(region, 1);
    }

    public final String refBases() { return mRefBases; }
    public void setRefBases(final String bases) { mRefBases = bases; }

    public double averageDepth()
    {
        int depthTotal = Arrays.stream(mRefBasesMatched).sum();
        return depthTotal/(double)mRefBasesMatched.length;
    }

    public String toString()
    {
        return String.format("%s: %s:%d -> %d", Id, chromosome(), start(), end());
    }
}
