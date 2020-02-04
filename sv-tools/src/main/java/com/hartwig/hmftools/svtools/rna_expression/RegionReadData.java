package com.hartwig.hmftools.svtools.rna_expression;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

public class RegionReadData
{
    public GenomeRegion Region;

    public final String Id;
    private final List<String> mRefRegions; // identifiers for this region, eg transcriptId

    private String mRefBases;
    private int[] mRefBasesMatched;
    private int[] mMatchTypeCounts;
    private final Map<RegionReadData,Integer> mLinkedRegions; // count of reads covering this region and another next to it
    private List<RegionReadData> mPreRegions; // references to adjacent regions with a lower position
    private List<RegionReadData> mPostRegions;

    public static final int MATCH_TYPE_NONE = -1;
    public static final int MATCH_TYPE_EXON_BOUNDARY = 0; // reads stopping at an exon boundary
    public static final int MATCH_TYPE_WITHIN_EXON = 1; // read fully contained within the exon
    public static final int MATCH_TYPE_SPAN_EXON_BOUNDARY = 2; //
    public static final int MATCH_TYPE_SPLICE_JUNCTION = 3; // read correctly ends at exon boundary of 2 adjacent exons
    public static final int MATCH_TYPE_UNSPLICED = 4; // reads spanning to unmapped regions where adjacent regions exist
    public static final int MATCH_TYPE_EXON_CHIMERIC = 5;
    public static final int MATCH_TYPE_INTRONIC = 6;
    public static final int MATCH_TYPE_MAX = MATCH_TYPE_INTRONIC + 1;

    public RegionReadData(final GenomeRegion region, final String id)
    {
        Region = region;
        Id = id;

        mRefRegions = Lists.newArrayList();

        mRefBases = "";
        mMatchTypeCounts = new int[MATCH_TYPE_MAX];

        mPreRegions = Lists.newArrayList();
        mPostRegions = Lists.newArrayList();
        mLinkedRegions = Maps.newHashMap();
    }

    public String chromosome() { return Region.chromosome(); }
    public long start() { return Region.start(); }
    public long end() { return Region.end(); }

    public void resetRegionBounds(long posStart, long posEnd)
    {
        Region = GenomeRegions.create(Region.chromosome(), posStart, posEnd);
    }

    public final List<String> getRefRegions() { return mRefRegions; }
    public void addRefRegion(final String ref)
    {
        if(!mRefRegions.contains(ref))
            mRefRegions.add(ref);
    }

    public void addMatchedRead(int matchType) { ++mMatchTypeCounts[matchType]; }
    public int matchedReadCount(int matchType) { return mMatchTypeCounts[matchType]; }

    public final String refBases() { return mRefBases; }

    public void setRefBases(final String bases)
    {
        mRefBases = bases;
        mRefBasesMatched = new int[(int)mRefBases.length()+1];
    }

    public int length() { return mRefBases.length(); }
    public int[] refBasesMatched() { return mRefBasesMatched; }

    public List<RegionReadData> getPreRegions() { return mPreRegions; }
    public List<RegionReadData> getPostRegions() { return mPostRegions; }

    public void addPreRegion(final RegionReadData region)
    {
        if(!mPreRegions.contains(region))
            mPreRegions.add(region);
    }

    public void addPostRegion(final RegionReadData region)
    {
        if(!mPostRegions.contains(region))
            mPostRegions.add(region);
    }

    public final Map<RegionReadData, Integer> getLinkedRegions() { return mLinkedRegions; }

    public void addLinkedRegion(final RegionReadData region)
    {
        Integer linkCount = mLinkedRegions.get(region);

        if(linkCount != null)
            mLinkedRegions.put(region, linkCount + 1);
        else
            mLinkedRegions.put(region, 1);
    }

    public double averageDepth()
    {
        if(mRefBasesMatched == null)
            return 0;

        int depthTotal = Arrays.stream(mRefBasesMatched).sum();
        return depthTotal/(double)mRefBasesMatched.length;
    }

    public int baseCoverage(int minReadCount)
    {
        if(mRefBasesMatched == null)
            return 0;

        // percent of bases covered by a read
        return (int)Arrays.stream(mRefBasesMatched).filter(x -> x >= minReadCount).count();
    }

    public String toString()
    {
        return String.format("%s: %s:%d -> %d refs(%d) reads(%d)", Id, chromosome(), start(), end(), mRefRegions.size(), mMatchTypeCounts);
    }

    public void clearState()
    {
        if(mRefBasesMatched != null)
        {
            for (int i = 0; i < mRefBasesMatched.length; ++i)
                mRefBasesMatched[i] = 0;
        }

        mLinkedRegions.clear();

        for(int i = 0; i < mMatchTypeCounts.length; ++i)
            mMatchTypeCounts[i] = 0;
    }
}
