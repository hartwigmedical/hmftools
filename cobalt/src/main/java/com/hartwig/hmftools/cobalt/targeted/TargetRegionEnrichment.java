package com.hartwig.hmftools.cobalt.targeted;

import java.util.Objects;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class TargetRegionEnrichment
{
    public final Chromosome mChromosome;
    public final int Position;
    public final double Enrichment;

    public TargetRegionEnrichment(final Chromosome mChromosome, final int position, final double enrichment)
    {
        this.mChromosome = mChromosome;
        Position = position;
        Enrichment = enrichment;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final TargetRegionEnrichment that = (TargetRegionEnrichment) o;
        return Position == that.Position && Objects.equals(mChromosome, that.mChromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mChromosome, Position);
    }

    @Override
    public String toString()
    {
        return "TargetRegionEnrichment{" +
                "mChromosome=" + mChromosome +
                ", Position=" + Position +
                ", Enrichment=" + Enrichment +
                '}';
    }
}
