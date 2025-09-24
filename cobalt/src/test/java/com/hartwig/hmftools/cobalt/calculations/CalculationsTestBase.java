package com.hartwig.hmftools.cobalt.calculations;

import static org.immutables.value.internal.$guava$.collect.$ImmutableList.of;

import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

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

    ChrBaseRegion cbr(Chromosome chromosome, int start, int end)
    {
        return new ChrBaseRegion(chromosome.contig(), start, end);
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
