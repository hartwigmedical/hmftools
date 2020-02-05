package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.MATCH_TYPE_INTRONIC;

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

    // summary results
    private int mTotalReadCount;
    private final Map<String,Integer> mTranscriptReadCounts;
    private final List<TranscriptResults> mTranscriptResults;

    public GeneReadData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;

        mExonRegions = Lists.newArrayList();
        mIntronRegions = Lists.newArrayList();
        mFragmentLengths = Lists.newArrayList();

        mTranscriptResults = Lists.newArrayList();
        mTotalReadCount = 0;
        mTranscriptReadCounts = Maps.newHashMap();
    }

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

        return intronRegions.stream().mapToInt(x -> x.matchedReadCount(MATCH_TYPE_INTRONIC)).sum();
    }

    public int totalReadCount() { return mTotalReadCount; }
    public void addReadCount() { ++mTotalReadCount; }
    public void setTotalReadCount(int count) { mTotalReadCount = count; }

    public final Map<String,Integer> getTranscriptReadCounts() { return mTranscriptReadCounts; }
    public int getTranscriptReadCount(final String trans)
    {
        Integer count = mTranscriptReadCounts.get(trans);
        return count != null ? count : 0;
    }

    public void addTranscriptReadMatch(final String trans)
    {
        Integer count = mTranscriptReadCounts.get(trans);
        if(count == null)
            mTranscriptReadCounts.put(trans, 1);
        else
            mTranscriptReadCounts.put(trans, count + 1);
    }

    public List<Integer> getFragmentLengths() { return mFragmentLengths; }
    public void addFragmentLength(int length) { mFragmentLengths.add(length); }

    public List<TranscriptResults> getTranscriptResults() { return mTranscriptResults; }
}
