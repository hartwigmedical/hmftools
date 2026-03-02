package com.hartwig.hmftools.amber.contamination;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.amber.TumorBAF;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;

public class PositionBasedFrequencySupplier implements TumorOnlyContaminationAnalysis.VariantFrequencySupplier
{
    private final Map<GenomePositionImpl, Double> FrequencyMap = new HashMap<>();

    public PositionBasedFrequencySupplier(Collection<AmberSite> sites)
    {
        for(AmberSite baf : sites)
        {
            FrequencyMap.put(new GenomePositionImpl(baf.chromosome(), baf.position()), baf.VariantAlleleFrequency);
        }
    }

    @Override
    public double minorAlleleFrequencyInPopulation(final TumorBAF tumorBAF)
    {
        return FrequencyMap.getOrDefault(new GenomePositionImpl(tumorBAF.chromosome(), tumorBAF.position()), 0.5);
    }
}
