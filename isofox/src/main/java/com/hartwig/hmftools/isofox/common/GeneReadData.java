package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.isofox.common.RegionReadData.generateExonicRegions;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.isofox.common.GeneMatchType.typeAsInt;
import static com.hartwig.hmftools.isofox.common.RnaUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsOverlap;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GeneReadData
{
    public final EnsemblGeneData GeneData;

    private final List<RegionReadData> mExonRegions; // set of unique exons ie with differing start and end positions
    private final List<long[]> mCommonExonicRegions; // merge any overlapping exons, to form a set of exonic regions for the gene

    private final List<TranscriptData> mTranscripts;
    private final long[] mTranscriptsRange;

    // summary results
    private final Map<Integer,int[][]> mTranscriptReadCounts; // count of fragments support types for each transcript, and whether unique
    private final Map<String,Double> mTranscriptAllocations; // results from the expected rate vs counts fit routine
    private double mFitResiduals;

    private final int[] mFragmentCounts;

    private static final Logger LOGGER = LogManager.getLogger(GeneReadData.class);

    public GeneReadData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;

        mExonRegions = Lists.newArrayList();
        mTranscripts = Lists.newArrayList();
        mCommonExonicRegions = Lists.newArrayList();

        mFragmentCounts = new int[typeAsInt(GeneMatchType.MAX)];
        mTranscriptReadCounts = Maps.newHashMap();
        mTranscriptAllocations = Maps.newHashMap();
        mFitResiduals = 0;

        mTranscriptsRange = new long[SE_PAIR];
    }

    public String name() { return GeneData.GeneName;}

    public final List<TranscriptData> getTranscripts() { return mTranscripts; }
    public void setTranscripts(final List<TranscriptData> transDataList) { mTranscripts.addAll(transDataList); }

    public final List<RegionReadData> getExonRegions() { return mExonRegions; }

    public final long[] getTranscriptsRange() { return mTranscriptsRange; }

    public void generateRegions()
    {
        // FIXME: still required?
        generateExonicRegions(GeneData.Chromosome, mExonRegions, mTranscripts);

        for(final TranscriptData transData : mTranscripts)
        {
            mTranscriptsRange[SE_END] = max(transData.TransEnd, mTranscriptsRange[SE_END]);

            if (mTranscriptsRange[SE_START] == 0 || transData.TransStart < mTranscriptsRange[SE_START])
                mTranscriptsRange[SE_START] = transData.TransStart;
        }

        generateCommonExonicRegions(mExonRegions, mCommonExonicRegions);
    }

    public final int[] getCounts() { return mFragmentCounts; }
    public void addCount(GeneMatchType type, int count) { mFragmentCounts[typeAsInt(type)] += count; }

    public static final int TRANS_COUNT = 0;
    public static final int UNIQUE_TRANS_COUNT = 1;

    public int[][] getTranscriptReadCount(final String trans)
    {
        int[][] counts = mTranscriptReadCounts.get(trans);
        return counts != null ? counts : new int[FragmentMatchType.MAX_FRAG_TYPE][UNIQUE_TRANS_COUNT+1];
    }

    public void addTranscriptReadMatch(int transId, boolean isUnique, FragmentMatchType type)
    {
        int[][] counts = mTranscriptReadCounts.get(transId);
        if(counts == null)
        {
            counts = new int[FragmentMatchType.MAX_FRAG_TYPE][UNIQUE_TRANS_COUNT+1];
            mTranscriptReadCounts.put(transId,  counts);
        }

        if(isUnique)
        {
            ++counts[FragmentMatchType.typeAsInt(type)][UNIQUE_TRANS_COUNT];
        }

        ++counts[FragmentMatchType.typeAsInt(type)][TRANS_COUNT];
    }

    public Map<String,Double> getTranscriptAllocations() { return mTranscriptAllocations; }

    public void setFitResiduals(double residuals) { mFitResiduals = residuals; }
    public double getFitResiduals() { return mFitResiduals; }

    public double getTranscriptAllocation(final String transName)
    {
        Double allocation = mTranscriptAllocations.get(transName);
        return allocation != null ? allocation : 0;
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

    public List<long[]> getCommonExonicRegions() { return mCommonExonicRegions; }

    public long calcExonicRegionLength()
    {
        return mCommonExonicRegions.stream().mapToLong(x -> x[SE_END] - x[SE_START]).sum();
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

    @VisibleForTesting
    public void clearCounts()
    {
        for(int i = 0; i < mFragmentCounts.length; ++i)
            mFragmentCounts[i] = 0;
    }

}
