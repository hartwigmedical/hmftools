package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.INTRONIC;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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
    private final List<Integer> mFragmentLengths;

    private final List<TranscriptData> mTranscripts;

    // summary results
    private final Map<String,int[]> mTranscriptReadCounts; // count of fragments supporting the transcript, and whether unique
    private final List<TranscriptResults> mTranscriptResults;

    private final int[] mFragmentCounts;

    public static final int GC_TOTAL = 0;
    public static final int GC_TRANS_SUPPORTING = 1;
    public static final int GC_ALT = 2;
    public static final int GC_INTRONIC = 3;
    public static final int GC_READ_THROUGH = 4;

    public GeneReadData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;

        mExonRegions = Lists.newArrayList();
        mIntronRegions = Lists.newArrayList();
        mFragmentLengths = Lists.newArrayList();
        mTranscripts = Lists.newArrayList();

        mTranscriptResults = Lists.newArrayList();
        mFragmentCounts = new int[GC_READ_THROUGH+1];
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

    public int getIntronRegionCounts(final TranscriptData transData, boolean restrictToUnique)
    {
        final List<RegionReadData> intronRegions = mIntronRegions.stream()
                .filter(x -> x.getRefRegions().contains(transData.TransName))
                .filter(x -> !restrictToUnique || x.getRefRegions().size() == 1)
                .collect(Collectors.toList());

        return intronRegions.stream().mapToInt(x -> x.matchedReadCount(INTRONIC)).sum();
    }

    public final int[] getCounts() { return mFragmentCounts; }
    public void addCount(int type, int count) { mFragmentCounts[type] += count; }

    public static final int TRANS_COUNT = 0;
    public static final int UNIQUE_TRANS_COUNT = 1;

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

    public List<Integer> getFragmentLengths() { return mFragmentLengths; }
    public void addFragmentLength(int length) { mFragmentLengths.add(length); }

    public List<TranscriptResults> getTranscriptResults() { return mTranscriptResults; }
}
