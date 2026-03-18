package com.hartwig.hmftools.amber.purity;

import java.util.Collection;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

public class GnomadFrequencyAnalysis
{
    private final GnomadFrequencySupplier GnomadData;

    public GnomadFrequencyAnalysis(final GnomadFrequencySupplier gnomadFrequencySupplier)
    {
        GnomadData = gnomadFrequencySupplier;
    }

    public double getMeanFrequency(final Collection<? extends GenomePosition> positions)
    {
        return positions.stream().mapToDouble(pos -> GnomadData.getFrequency(pos.chromosome(), pos.position())).average().orElse(0);
    }
}
