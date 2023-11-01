package com.hartwig.hmftools.sage.coverage;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

public class GeneCoverage
{
    private final String mChromosome;
    private final String mGene;
    private final List<ExonCoverage> mExonCoverage;
    private final int mMinPosition;
    private final int mMaxPosition;
    private final int mMinMapQuality;

    public GeneCoverage(final String geneName, final List<ExonCoverage> exons, int minMapQuality)
    {
        mChromosome = exons.get(0).chromosome();
        mGene = geneName;
        mExonCoverage = exons;
        mMinMapQuality = minMapQuality;

        int tmpMin = exons.get(0).start();
        int tmpMax = exons.get(0).end();

        for(ExonCoverage exon : exons)
        {
            tmpMin = Math.min(tmpMin, exon.start());
            tmpMax = Math.max(tmpMax, exon.end());
        }

        mMinPosition = tmpMin;
        mMaxPosition = tmpMax;
    }

    public List<ExonCoverage> exonCoverage() { return mExonCoverage; }
    public String geneName() { return mGene; }
    public int minPosition() { return mMinPosition; }
    public int maxPosition() { return mMaxPosition; }

    public GeneCoverage clone()
    {
        List<ExonCoverage> exonCoverages = Lists.newArrayList();
        GeneCoverage coverage = new GeneCoverage(mGene, exonCoverages, mMinMapQuality);
        mExonCoverage.stream().forEach(x -> exonCoverages.add(new ExonCoverage(x.region(), x.exonRank())));
        return coverage;
    }

    public void processRead(final SAMRecord record)
    {
        if(record.getMappingQuality() < mMinMapQuality)
            return;

        if(!record.getContig().equals(mChromosome))
            return;

        int readStartPos = record.getAlignmentStart();
        int readEndPos = record.getAlignmentEnd();

        if(!positionsOverlap(readStartPos, readEndPos, mMinPosition, mMaxPosition))
            return;

        mExonCoverage.forEach(x -> x.processRead(readStartPos, readEndPos));
    }

    public String chromosome()
    {
        return mChromosome;
    }
}
