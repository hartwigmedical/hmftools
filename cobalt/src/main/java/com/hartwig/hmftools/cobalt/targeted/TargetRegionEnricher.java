package com.hartwig.hmftools.cobalt.targeted;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.calculations.CobaltCalculation;
import com.hartwig.hmftools.cobalt.calculations.CobaltWindow;
import com.hartwig.hmftools.cobalt.count.ReadDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class TargetRegionEnricher implements CobaltCalculation.TargetRegions
{
    private final ListMultimap<Chromosome, TargetRegionEnrichment> mEnrichments;

    public TargetRegionEnricher(final ListMultimap<Chromosome, TargetRegionEnrichment> mEnrichments)
    {
        this.mEnrichments = mEnrichments;
    }

    @Override
    public boolean isInTargetRegions(final Chromosome chromosome, final CobaltWindow window)
    {
        // todo do something more efficient
        return mEnrichments.get(chromosome)
                .stream()
                .anyMatch(targetRegionEnrichment -> targetRegionEnrichment.Position == window.Position);

    }

    @Override
    public boolean applyFinalNormalisation()
    {
        return true;
    }

    @Override
    public double enrichmentQuotient(final Chromosome chromosome, final ReadDepth readDepth)
    {
        // todo do something more efficient
        return mEnrichments.get(chromosome)
                .stream()
                .filter(targetRegionEnrichment -> targetRegionEnrichment.Position == readDepth.StartPosition)
                .findFirst()
                .orElseThrow().Enrichment;
    }
}
