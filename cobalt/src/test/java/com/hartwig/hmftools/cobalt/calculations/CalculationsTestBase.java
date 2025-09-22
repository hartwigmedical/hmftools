package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.immutables.value.internal.$guava$.collect.$ImmutableList.of;

import java.util.List;

import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCProfile;

import org.junit.Assert;
import org.junit.Test;

public class CalculationsTestBase
{
    BamRatio br(Chromosome chromosome, int pos, double depth, double gc, boolean onTarget)
    {
        return new BamRatio(chromosome, dr(chromosome, pos, depth, gc), onTarget);
    }

    DepthReading dr(Chromosome chromosome, int position, double depth, double gc)
    {
        return new DepthReading(chromosome.contig(), position, depth, gc);
    }

    GCProfile gcProfile(Chromosome chromosome, int position, double gc, double mappablePercentage)
    {
        return ImmutableGCProfile.builder()
                .chromosome(chromosome.contig())
                .start(position)
                .end(position + 1000)
                .nonNPercentage(100.0)
                .mappablePercentage(mappablePercentage)
                .gcContent(gc)
                .build();
    }
}
