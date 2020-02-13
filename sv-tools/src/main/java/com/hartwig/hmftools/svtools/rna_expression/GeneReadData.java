package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.MAX;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.typeAsInt;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.INTRONIC;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

public class GeneReadData
{
    public final EnsemblGeneData GeneData;

    private final List<RegionReadData> mExonRegions;
    private final List<RegionReadData> mIntronRegions;

    private final List<TranscriptData> mTranscripts;

    // summary results
    private final Map<String,int[][]> mTranscriptReadCounts; // count of fragments support types for each transcript, and whether unique
    private final List<TranscriptResults> mTranscriptResults;

    public static final int TC_SPLICE = 0;
    public static final int TC_SHORT = 1;
    public static final int TC_LONG = 2;

    private final int[] mFragmentCounts;

    public GeneReadData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;

        mExonRegions = Lists.newArrayList();
        mIntronRegions = Lists.newArrayList();
        mTranscripts = Lists.newArrayList();

        mTranscriptResults = Lists.newArrayList();
        mFragmentCounts = new int[typeAsInt(GeneMatchType.MAX)];
        mTranscriptReadCounts = Maps.newHashMap();
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

    public final List<RegionReadData> getIntronRegions() { return mIntronRegions; }

    public RegionReadData createOrFindIntronRegion(long intronStart, long intronEnd)
    {
        RegionReadData intronRegion = mIntronRegions.stream()
                .filter(x -> x.Region.start() <= intronStart && intronEnd <= x.Region.end())
                .findFirst().orElse(null);

        if(intronRegion == null)
        {
            GenomeRegion region = GenomeRegions.create(GeneData.Chromosome, intronStart, intronEnd);
            intronRegion = new RegionReadData(region);
            mIntronRegions.add(intronRegion);
        }
        else
        {
            // tighten intron region bounds
            if(intronStart > intronRegion.start() || intronEnd < intronRegion.end())
            {
                intronRegion.resetRegionBounds(max(intronStart, intronRegion.start()), min(intronEnd, intronRegion.end()));
            }
        }

        return intronRegion;
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

    public List<TranscriptResults> getTranscriptResults() { return mTranscriptResults; }

    @VisibleForTesting
    public void clearCounts()
    {
        for(int i = 0; i < mFragmentCounts.length; ++i)
            mFragmentCounts[i] = 0;
    }

}
