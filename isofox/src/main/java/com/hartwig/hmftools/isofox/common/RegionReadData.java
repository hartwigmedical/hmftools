package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.isofox.common.GeneReadData.TRANS_COUNT;
import static com.hartwig.hmftools.isofox.common.GeneReadData.UNIQUE_TRANS_COUNT;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

public class RegionReadData implements Comparable< RegionReadData>
{
    public GenomeRegion Region;

    private final List<String> mRefRegions; // identifiers for this region, eg transcriptId

    private String mRefBases;
    private int[] mRefBasesMatched;
    private boolean[] mUniqueRefBases; // bases not covered by any other region

    private final Map<Integer,int[]> mTranscriptReadCounts; // count of reads which support this region and a specific transcript
    private final Map<Integer,int[][]> mTranscriptJunctionCounts; // count of reads which support each exon junction and a specific transcript

    private List<RegionReadData> mPreRegions; // references to adjacent regions with a lower position
    private List<RegionReadData> mPostRegions;

    public RegionReadData(final GenomeRegion region)
    {
        Region = region;

        mRefRegions = Lists.newArrayList();

        mRefBases = "";

        mPreRegions = Lists.newArrayList();
        mPostRegions = Lists.newArrayList();

        mTranscriptReadCounts = Maps.newHashMap();
        mTranscriptJunctionCounts = Maps.newHashMap();
    }

    public String chromosome() { return Region.chromosome(); }
    public long start() { return Region.start(); }
    public long end() { return Region.end(); }

    public static final int NO_EXON = -1;
    public static final int TRANS_ID = 0;
    public static final int TRANS_NAME = 1;
    public static final int EXON_RANK = 2;
    public static final String TRANS_REF_DELIM = ":";

    private static String formExonRefId(int transId, final String transName, int exonRank)
    {
        return String.format("%d:%s:%d", transId, transName, exonRank);
    }

    public static int extractTransId(final String ref) { return Integer.valueOf(ref.split(TRANS_REF_DELIM)[TRANS_ID]); }
    public static String extractTransName(final String ref) { return ref.split(TRANS_REF_DELIM)[TRANS_NAME]; }
    public static int extractExonRank(final String ref) { return Integer.valueOf(ref.split(TRANS_REF_DELIM)[EXON_RANK]); }

    public int getExonRank(final int transId)
    {
        for(String transRef : mRefRegions)
        {
            if(extractTransId(transRef) == transId)
                return extractExonRank(transRef);
        }

        return NO_EXON;
    }

    public final List<String> getRefRegions() { return mRefRegions; }

    public boolean hasTransId(final int transId)
    {
        return mRefRegions.stream().anyMatch(x -> extractTransId(x) == transId);
    }

    public void addExonRef(int transId, final String transName, int exonRank)
    {
        if(!mRefRegions.contains(transId))
            mRefRegions.add(formExonRefId(transId, transName, exonRank));
    }

    public final String refBases() { return mRefBases; }

    public void setRefBases(final String bases)
    {
        mRefBases = bases;
        mRefBasesMatched = new int[mRefBases.length()];
        mUniqueRefBases = new boolean[mRefBases.length()];

        for(int i = 0; i < mUniqueRefBases.length; ++i)
            mUniqueRefBases[i] = true;
    }

    public int length() { return mRefBases.length(); }
    public int[] refBasesMatched() { return mRefBasesMatched; }
    public boolean[] uniqueRefBases() { return mUniqueRefBases; }

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

    public int[] getTranscriptReadCount(int transId)
    {
        int[] counts = mTranscriptReadCounts.get(transId);
        return counts != null ? counts : new int[]{0, 0};
    }

    public void addTranscriptReadMatch(int transId, boolean isUnique)
    {
        int[] counts = mTranscriptReadCounts.get(transId);
        if(counts == null)
        {
            counts = new int[2];
            mTranscriptReadCounts.put(transId,  counts);
        }

        if(isUnique)
            ++counts[UNIQUE_TRANS_COUNT];

        ++counts[TRANS_COUNT];
    }

    public int[] getTranscriptJunctionMatchCount(int transId, int seIndex)
    {
        int[][] counts = mTranscriptJunctionCounts.get(transId);
        return counts != null ? counts[seIndex] : new int[]{0, 0};
    }

    public void addTranscriptJunctionMatch(int transId, int seIndex, boolean isUnique)
    {
        int[][] counts = mTranscriptJunctionCounts.get(transId);
        if(counts == null)
        {
            counts = new int[SE_PAIR][2];
            mTranscriptJunctionCounts.put(transId, counts);
        }

        if(isUnique)
            ++counts[seIndex][UNIQUE_TRANS_COUNT];

        ++counts[seIndex][TRANS_COUNT];
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

        return (int)Arrays.stream(mRefBasesMatched).filter(x -> x >= minReadCount).count();
    }

    public static void findUniqueBases(final List<RegionReadData> regions)
    {
        for(int i = 0; i < regions.size() - 1; ++i)
        {
            final RegionReadData region1 = regions.get(i);

            for(int j = i + 1; j < regions.size(); ++j)
            {
                final RegionReadData region2 = regions.get(j);
                region1.markNonUniqueBases(region2);
            }
        }
    }

    public void markNonUniqueBases(final RegionReadData other)
    {
        if(start() > other.end() || end() < other.start())
            return;

        int overlapBases = (int)(min(end(), other.end()) - max(start(), other.start())) + 1;

        int thisOffset = (int)max(other.start() - start(), 0);
        int otherOffset = (int)max(start() - other.start(), 0);

        for(int i = 0; i < overlapBases; ++i)
        {
            mUniqueRefBases[i + thisOffset] = false;
            other.uniqueRefBases()[i + otherOffset] = false;
        }
    }

    public int uniqueBaseCount()
    {
        int count = 0;
        for(int i = 0; i < mUniqueRefBases.length; ++i)
        {
            if(mUniqueRefBases[i])
                ++count;
        }

        return count;
    }

    public int uniqueBaseCoverage(int minReadCount)
    {
        if(mRefBasesMatched == null || mUniqueRefBases == null)
            return 0;

        int count = 0;
        for(int i = 0; i < mRefBasesMatched.length; ++i)
        {
            if(mUniqueRefBases[i] && mRefBasesMatched[i] >= minReadCount)
                ++count;
        }

        return count;
    }

    public int uniqueBaseTotalDepth()
    {
        if(mRefBasesMatched == null || mUniqueRefBases == null)
            return 0;

        int total = 0;
        for(int i = 0; i < mRefBasesMatched.length; ++i)
        {
            if(mUniqueRefBases[i])
                total += mRefBasesMatched[i];
        }

        return total;
    }

    public String toString()
    {
        int sjReads = mTranscriptJunctionCounts.values().stream()
                .mapToInt(x -> x[SE_START][TRANS_COUNT] + x[SE_END][TRANS_COUNT]).sum();

        int reads = mTranscriptReadCounts.values().stream().mapToInt(x -> x[TRANS_COUNT]).sum();

        return String.format("%s %s:%d -> %d refs(%d) %s",
                !mRefRegions.isEmpty() ? mRefRegions.get(0) : "unknown", chromosome(), start(), end(), mRefRegions.size(),
                mRefBases != null ? String.format("reads(%d sj=%d)",reads, sjReads) : "intron");
    }

    public void clearState()
    {
        if(mRefBasesMatched != null)
        {
            for (int i = 0; i < mRefBasesMatched.length; ++i)
                mRefBasesMatched[i] = 0;
        }
    }

    @Override
    public int compareTo(RegionReadData other)
    {
        return (int)(start() - other.start());
    }

}
