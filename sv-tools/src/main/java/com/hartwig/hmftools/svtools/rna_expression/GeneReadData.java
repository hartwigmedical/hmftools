package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.typeAsInt;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.deriveCommonRegions;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

public class GeneReadData
{
    public final EnsemblGeneData GeneData;

    private final List<RegionReadData> mExonRegions;
    private final List<long[]> mCommonExonicRegions;

    private final List<TranscriptData> mTranscripts;
    private final long[] mTranscriptsRange;

    // summary results
    private final Map<String, int[][]> mTranscriptReadCounts; // count of fragments support types for each transcript, and whether unique
    private final Map<String, Double> mTranscriptAllocations;
    private final List<TranscriptResults> mTranscriptResults;

    public static final int TC_SPLICED = 0;
    public static final int TC_SHORT = 1;
    public static final int TC_LONG = 2;
    public static final int TC_UNSPLICED = 3;
    public static final int TC_MAX = 4;

    private final int[] mFragmentCounts;

    public GeneReadData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;

        mExonRegions = Lists.newArrayList();
        mTranscripts = Lists.newArrayList();
        mCommonExonicRegions = Lists.newArrayList();

        mTranscriptResults = Lists.newArrayList();
        mFragmentCounts = new int[typeAsInt(GeneMatchType.MAX)];
        mTranscriptReadCounts = Maps.newHashMap();
        mTranscriptAllocations = Maps.newHashMap();

        mTranscriptsRange = new long[SE_PAIR];
    }

    public final List<TranscriptData> getTranscripts() { return mTranscripts; }
    public void setTranscripts(final List<TranscriptData> transDataList) { mTranscripts.addAll(transDataList); }

    public final List<RegionReadData> getExonRegions() { return mExonRegions; }

    public void addExonRegion(final RegionReadData exonRegionData)
    {
        if(hasExonRegion(exonRegionData.Region.start(), exonRegionData.Region.end()))
            return;

        mExonRegions.add(exonRegionData);
    }

    public boolean hasExonRegion(long posStart, long posEnd)
    {
        return mExonRegions.stream().anyMatch(x -> x.Region.start() == posStart && x.Region.end() == posEnd);
    }

    public RegionReadData findExonRegion(long posStart, long posEnd)
    {
        return mExonRegions.stream()
                .filter(x -> x.Region.start() == posStart && x.Region.end() == posEnd)
                .findFirst().orElse(null);
    }

    public final long[] getTranscriptsRange() { return mTranscriptsRange; }

    public void generateExonicRegions()
    {
        // form a genomic region for each unique exon amongst the transcripts
        for(final TranscriptData transData : mTranscripts)
        {
            RegionReadData prevRegionReadData = null;

            for(int i = 0; i < transData.exons().size(); ++ i)
            {
                ExonData exon = transData.exons().get(i);

                RegionReadData exonReadData = findExonRegion(exon.ExonStart, exon.ExonEnd);

                if (exonReadData == null)
                {
                    GenomeRegion region = GenomeRegions.create(GeneData.Chromosome, exon.ExonStart, exon.ExonEnd);
                    exonReadData = new RegionReadData(region);
                    addExonRegion(exonReadData);
                }

                exonReadData.addExonRef(transData.TransName, exon.ExonRank);

                if(prevRegionReadData != null)
                {
                    prevRegionReadData.addPostRegion(exonReadData);
                    exonReadData.addPreRegion(prevRegionReadData);
                }

                prevRegionReadData = exonReadData;
            }

            mTranscriptsRange[SE_END] = max(transData.TransEnd, mTranscriptsRange[SE_END]);

            if(mTranscriptsRange[SE_START] == 0 || transData.TransStart < mTranscriptsRange[SE_START])
                mTranscriptsRange[SE_START] = transData.TransStart;
        }

        generateCommonExonicRegions();
    }

    public final int[] getCounts() { return mFragmentCounts; }
    public void addCount(GeneMatchType type, int count) { mFragmentCounts[typeAsInt(type)] += count; }

    public static final int TRANS_COUNT = 0;
    public static final int UNIQUE_TRANS_COUNT = 1;

    public int[][] getTranscriptReadCount(final String trans)
    {
        int[][] counts = mTranscriptReadCounts.get(trans);
        return counts != null ? counts : new int[TC_LONG+1][UNIQUE_TRANS_COUNT+1];
    }

    public void addTranscriptReadMatch(final String trans, boolean isUnique, int type)
    {
        int[][] counts = mTranscriptReadCounts.get(trans);
        if(counts == null)
        {
            counts = new int[TC_LONG+1][UNIQUE_TRANS_COUNT+1];
            mTranscriptReadCounts.put(trans,  counts);
        }

        if(isUnique)
        {
            ++counts[type][UNIQUE_TRANS_COUNT];
        }

        ++counts[type][TRANS_COUNT];
    }

    public Map<String,Double> getTranscriptAllocations() { return mTranscriptAllocations; }

    public List<TranscriptResults> getTranscriptResults() { return mTranscriptResults; }

    private void generateCommonExonicRegions()
    {
        if(mExonRegions.isEmpty())
            return;

        List<long[]> commonExonicRegions = Lists.newArrayList(new long[] {mExonRegions.get(0).start(), mExonRegions.get(0).end()});

        for(int i = 1; i < mExonRegions.size(); ++i)
        {
            List<long[]> nextRegion = Lists.newArrayList(new long[] {mExonRegions.get(i).start(), mExonRegions.get(i).end()});
            commonExonicRegions = deriveCommonRegions(commonExonicRegions, nextRegion);
        }

        mCommonExonicRegions.addAll(commonExonicRegions);
    }

    public List<long[]> getCommonExonicRegions() { return mCommonExonicRegions; }

    public long calcExonicRegionLength()
    {
        return mCommonExonicRegions.stream().mapToLong(x -> x[SE_END] - x[SE_START]).sum();
    }

    public static String countsTypeToStr(int type)
    {
        switch(type)
        {
            case TC_SHORT: return "SHORT";
            case TC_LONG: return "LONG";
            case TC_UNSPLICED: return "UNSPLICED";
            case TC_SPLICED: return "SPLICED";
            default: return "UNKNOWN";
        }
    }

    @VisibleForTesting
    public void clearCounts()
    {
        for(int i = 0; i < mFragmentCounts.length; ++i)
            mFragmentCounts[i] = 0;
    }

}
