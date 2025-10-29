package com.hartwig.hmftools.common.bam.testutilities;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

public record ChromosomeWindow(HumanChromosome chromosome, int start, int end)
{
    public ChromosomeWindow next(int step)
    {
        return new ChromosomeWindow(chromosome, start + step, end + step);
    }

    public String chromosomeName()
    {
        return V38.versionedChromosome(chromosome);
    }

    public Pair<BaseRegion,BaseRegion> toACGTBasesRegion()
    {
        int length = end - start;
        int quadCount = length / 4;
        byte[] bases = "ACGT".repeat(Math.max(0, quadCount)).getBytes();
        return buildPairFromBases(bases);
    }

    public Pair<BaseRegion,BaseRegion> toRandomBasesRegion()
    {
        int length = end - start;
        byte[] bases = MockRefGenome.generateRandomBases(length).getBytes();
        return buildPairFromBases(bases);
    }

    public Pair<BaseRegion, BaseRegion> toBasesWithGivenGC(double gcRatio)
    {
        Preconditions.checkArgument(gcRatio >= 0.0 && gcRatio <= 1.0);
        int length = end - start;
        Preconditions.checkArgument(length > 0);
        Preconditions.checkArgument(length % 100 == 0);
        int gCount = (int) (length * gcRatio);
        String gPart = "G".repeat(gCount);
        String aPart = "A".repeat(length - gCount);
        byte[] bases = (gPart + aPart).getBytes();
        return buildPairFromBases(bases);
    }

    public Pair<BaseRegion,BaseRegion> toBaseRegionPair(RefGenomeSource refGenomeSource){
        return Pair.of(toBaseRegion(refGenomeSource), mateBaseRegion(refGenomeSource));
    }

    private BaseRegion toBaseRegion(RefGenomeSource refGenomeSource)
    {
        byte[] bases = refGenomeSource.getBases(chromosomeName(), start, end - 1);
        return new BaseRegion(chromosome, start, end, bases);
    }

    private BaseRegion mateBaseRegion(RefGenomeSource refGenomeSource)
    {
        int mateStart = end;
        int mateStop = mateStart + (end - start);
        byte[] mateBases = refGenomeSource.getBases(chromosomeName(), mateStart, mateStop - 1);
        return new BaseRegion(chromosome, mateStart, mateStop, mateBases);
    }

    @NotNull
    private ImmutablePair<BaseRegion, BaseRegion> buildPairFromBases(final byte[] bases)
    {
        BaseRegion left = new BaseRegion(chromosome, start, end, bases);
        int mateStart = end;
        int mateStop = mateStart + (end - start);
        BaseRegion right = new BaseRegion(chromosome, mateStart, mateStop, bases);
        return new ImmutablePair<>(left, right);
    }
}
