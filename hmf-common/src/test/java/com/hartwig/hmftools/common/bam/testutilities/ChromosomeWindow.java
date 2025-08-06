package com.hartwig.hmftools.common.bam.testutilities;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

public record ChromosomeWindow(int chromosome, int start, int end)
{
    public ChromosomeWindow next(int step)
    {
        return new ChromosomeWindow(chromosome, start + step, end + step);
    }

    public String chromosomeName()
    {
        return "chr" + (chromosome + 1);
    }

    public Pair<BaseRegion,BaseRegion> toRandomBasesRegion()
    {
        int length = end - start;
        byte[] bases = MockRefGenome.generateRandomBases(length).getBytes();
        BaseRegion left = new BaseRegion(chromosome, start, end, bases);
        int mateStart = end;
        int mateStop = mateStart + (end - start);
        BaseRegion right = new BaseRegion(chromosome, mateStart, mateStop, bases);
        return new ImmutablePair<>(left, right);
    }

    public Pair<BaseRegion,BaseRegion> toBaseRegionPair(RefGenomeSource refGenomeSource){
        return Pair.of(toBaseRegion(refGenomeSource), mateBaseRegion(refGenomeSource));
    }

    public BaseRegion toBaseRegion(RefGenomeSource refGenomeSource)
    {
        byte[] bases = refGenomeSource.getBases(chromosomeName(), start, end - 1);
        return new BaseRegion(chromosome, start, end, bases);
    }

    public BaseRegion mateBaseRegion(RefGenomeSource refGenomeSource)
    {
        int mateStart = end;
        int mateStop = mateStart + (end - start);
        byte[] mateBases = refGenomeSource.getBases(chromosomeName(), mateStart, mateStop - 1);
        return new BaseRegion(chromosome, mateStart, mateStop, mateBases);
    }
}
