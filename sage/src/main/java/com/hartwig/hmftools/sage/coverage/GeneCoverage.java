package com.hartwig.hmftools.sage.coverage;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.bed.NamedBed;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.math3.distribution.PoissonDistribution;

public class GeneCoverage
{
    private final String mChromosome;
    private final String mGene;
    private final List<ExonCoverage> mExonCoverage;
    private final int mMinPosition;
    private final int mMaxPosition;

    public GeneCoverage(final String geneName, final List<ExonCoverage> exons)
    {
        mChromosome = exons.get(0).chromosome();
        mGene = geneName;
        mExonCoverage = exons; // .stream().map(ExonCoverage::new).collect(Collectors.toList());

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
        GeneCoverage coverage = new GeneCoverage(mGene, exonCoverages);
        mExonCoverage.stream().forEach(x -> exonCoverages.add(new ExonCoverage(x.region(), x.exonRank())));
        return coverage;
    }

    public void processRead(final String chromosome, int readStartPos, int readEndPos)
    {
        if(!chromosome.equals(mChromosome) || !positionsOverlap(readStartPos, readEndPos, mMinPosition, mMaxPosition))
            return;

        mExonCoverage.forEach(x -> x.processRead(readStartPos, readEndPos));
    }

    public String chromosome()
    {
        return mChromosome;
    }

}
