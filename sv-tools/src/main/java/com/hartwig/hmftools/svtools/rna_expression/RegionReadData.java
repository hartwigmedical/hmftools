package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TRANS_COUNT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.UNIQUE_TRANS_COUNT;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

public class RegionReadData implements Comparable< RegionReadData>
{
    public GenomeRegion Region;

    private final List<String> mRefRegions; // identifiers for this region, eg transcriptId

    private String mRefBases;
    private int[] mRefBasesMatched;

    private final Map<String, int[]> mTranscriptReadCounts; // count of reads which support this region and a specific transcript
    private final Map<String, int[][]> mTranscriptJunctionCounts; // count of reads which support each exon junction and a specific transcript

    private int[] mMatchTypeCounts;
    private List<RegionReadData> mPreRegions; // references to adjacent regions with a lower position
    private List<RegionReadData> mPostRegions;

    public RegionReadData(final GenomeRegion region)
    {
        Region = region;

        mRefRegions = Lists.newArrayList();

        mRefBases = "";
        mMatchTypeCounts = new int[RegionMatchType.values().length];

        mPreRegions = Lists.newArrayList();
        mPostRegions = Lists.newArrayList();

        mTranscriptReadCounts = Maps.newHashMap();
        mTranscriptJunctionCounts = Maps.newHashMap();
    }

    public String chromosome() { return Region.chromosome(); }
    public long start() { return Region.start(); }
    public long end() { return Region.end(); }

    public void resetRegionBounds(long posStart, long posEnd)
    {
        Region = GenomeRegions.create(Region.chromosome(), posStart, posEnd);
    }

    public static final int NO_EXON = -1;
    public static final int TRANS_ID = 0;
    public static final int EXON_RANK = 1;

    private static String formExonRefId(final String transId, int exonRank) {  return String.format("%s:%d", transId, exonRank); }
    public static String extractTransId(final String ref) { return ref.split(":")[TRANS_ID]; }
    public static int extractExonRank(final String ref) { return Integer.valueOf(ref.split(":")[EXON_RANK]); }

    public int getExonRank(final String transId)
    {
        final String exonRefId = mRefRegions.stream().filter(x -> x.contains(transId)).findFirst().orElse(null);
        return exonRefId != null ? extractExonRank(exonRefId) : NO_EXON;
    }

    public final List<String> getRefRegions() { return mRefRegions; }

    public boolean hasTransId(final String transId)
    {
        return mRefRegions.stream().anyMatch(x -> x.contains(transId));
    }

    public void addExonRef(final String transId, int exonRank)
    {
        if(!mRefRegions.contains(transId))
            mRefRegions.add(formExonRefId(transId, exonRank));
    }

    public void addMatchedRead(RegionMatchType matchType) { ++mMatchTypeCounts[matchType.ordinal()]; }
    public int matchedReadCount(RegionMatchType matchType) { return mMatchTypeCounts[matchType.ordinal()]; }

    public final String refBases() { return mRefBases; }

    public void setRefBases(final String bases)
    {
        mRefBases = bases;
        mRefBasesMatched = new int[(int)mRefBases.length()];
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

    public int[] getTranscriptReadCount(final String trans)
    {
        int[] counts = mTranscriptReadCounts.get(trans);
        return counts != null ? counts : new int[]{0, 0};
    }

    public void addTranscriptReadMatch(final String trans, boolean isUnique)
    {
        int[] counts = mTranscriptReadCounts.get(trans);
        if(counts == null)
        {
            counts = new int[2];
            mTranscriptReadCounts.put(trans,  counts);
        }

        if(isUnique)
            ++counts[UNIQUE_TRANS_COUNT];

        ++counts[TRANS_COUNT];
    }

    public int[] getTranscriptJunctionMatchCount(final String trans, int seIndex)
    {
        int[][] counts = mTranscriptJunctionCounts.get(trans);
        return counts != null ? counts[seIndex] : new int[]{0, 0};
    }

    public void addTranscriptJunctionMatch(final String trans, int seIndex, boolean isUnique)
    {
        int[][] counts = mTranscriptJunctionCounts.get(trans);
        if(counts == null)
        {
            counts = new int[SE_PAIR][2];
            mTranscriptJunctionCounts.put(trans, counts);
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

        for(int i = 0; i < mMatchTypeCounts.length; ++i)
            mMatchTypeCounts[i] = 0;
    }

    @Override
    public int compareTo(RegionReadData other)
    {
        return (int)(start() - other.start());
    }

}
