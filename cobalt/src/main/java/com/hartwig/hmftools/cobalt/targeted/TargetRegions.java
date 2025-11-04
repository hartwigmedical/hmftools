package com.hartwig.hmftools.cobalt.targeted;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.normalisers.NoOpReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.UnityNormaliser;
import com.hartwig.hmftools.cobalt.consolidation.NoOpConsolidator;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class TargetRegions implements CobaltScope
{
    public interface ChromosomeData
    {
        int length(Chromosome chromosome);
    }

    private final Map<Chromosome, ArrayList<TargetRegionEnrichment>> mEnrichments = new HashMap<>();

    public TargetRegions(ListMultimap<Chromosome, TargetRegionEnrichment> enrichments, ChromosomeData chromosomeData)
    {
        for(Chromosome chromosome : enrichments.keySet())
        {
            int length = chromosomeData.length(chromosome);
            int numberOfSlots = length / 1000;
            ArrayList<TargetRegionEnrichment> enrichmentsForChromosome = new ArrayList<>(numberOfSlots);
            int position = 1;
            Map<Integer, TargetRegionEnrichment> positionToSuppliedItem = new HashMap<>();
            for(TargetRegionEnrichment enrichment : enrichments.get(chromosome))
            {
                positionToSuppliedItem.put(enrichment.Position, enrichment);
            }

            for(int i = 0; i < numberOfSlots; i++)
            {
                TargetRegionEnrichment suppliedItem = positionToSuppliedItem.get(position);
                enrichmentsForChromosome.add(i, suppliedItem);
                position += 1000;
            }
            mEnrichments.put(chromosome, enrichmentsForChromosome);
        }
    }

    @Override
    public ReadDepthStatisticsNormaliser medianByMeanNormaliser()
    {
        return new NoOpReadDepthStatisticsNormaliser();
    }

    @Override
    public ResultsNormaliser finalNormaliser()
    {
        return new UnityNormaliser();
    }

    @Override
    public ResultsConsolidator resultsConsolidator(final double medianReadDepth)
    {
        return new NoOpConsolidator();
    }

    @Override
    public boolean onTarget(final Chromosome chromosome, final int position)
    {
        return getEnrichment(chromosome, position) != null;
    }

    @Override
    public double enrichmentQuotient(final Chromosome chromosome, final int position)
    {
        TargetRegionEnrichment enrichment = getEnrichment(chromosome, position);
        return enrichment == null ? -1.0 : enrichment.Enrichment;
    }

    private TargetRegionEnrichment getEnrichment(final Chromosome chromosome, final int position)
    {
        if(mEnrichments.containsKey(chromosome))
        {
            return mEnrichments.get(chromosome).get(position / 1000);
        }
        return null;
    }
}
