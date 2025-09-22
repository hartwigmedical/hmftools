package com.hartwig.hmftools.cobalt.targeted;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.calculations.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.calculations.UnityNormaliser;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class TargetRegionEnricher implements TargetRegions
{
    public interface ChromosomeData
    {
        int length(Chromosome chromosome);
    }

    private final Map<Chromosome, ArrayList<TargetRegionEnrichment>> mEnrichments = new HashMap<>();

    public TargetRegionEnricher(ListMultimap<Chromosome, TargetRegionEnrichment> enrichments, ChromosomeData chromosomeData)
    {
        for (Chromosome chromosome : enrichments.keySet())
        {
            int length = chromosomeData.length(chromosome);
            int numberOfSlots = length / 1000;
            ArrayList<TargetRegionEnrichment> enrichmentsForChromosome = new ArrayList<>(numberOfSlots);
            int position = 1;
            Map<Integer, TargetRegionEnrichment> positionToSuppliedItem = new HashMap<>();
            for (TargetRegionEnrichment enrichment : enrichments.get(chromosome))
            {
                positionToSuppliedItem.put(enrichment.Position, enrichment);
            }

            for (int i = 0; i < numberOfSlots; i++)
            {
                TargetRegionEnrichment suppliedItem = positionToSuppliedItem.get(position);
                enrichmentsForChromosome.add(i, suppliedItem);
                position += 1000;
            }
            mEnrichments.put(chromosome, enrichmentsForChromosome);
        }
    }

    @Override
    public boolean onTarget(final Chromosome chromosome, final int position)
    {
        return getEnrichment(chromosome, position) != null;
    }

    @Override
    public ResultsNormaliser createNormaliser()
    {
        return new UnityNormaliser();
    }

    @Override
    public double enrichmentQuotient(final Chromosome chromosome, final DepthReading readDepth)
    {
        TargetRegionEnrichment enrichment = getEnrichment(chromosome, readDepth.StartPosition);
        return enrichment == null ? -1.0 : enrichment.Enrichment;
    }

    private TargetRegionEnrichment getEnrichment(final Chromosome chromosome, final int position)
    {
        if (mEnrichments.containsKey(chromosome))
        {
            return mEnrichments.get(chromosome).get(position / 1000);
        }
        return null;
    }
}
