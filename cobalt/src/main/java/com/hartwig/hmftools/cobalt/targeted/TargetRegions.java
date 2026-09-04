package com.hartwig.hmftools.cobalt.targeted;

import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.WINDOW_SIZE;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.cobalt.normalisers.NoOpReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.ResultsNormaliser;
import com.hartwig.hmftools.cobalt.normalisers.UnityNormaliser;
import com.hartwig.hmftools.cobalt.consolidation.NoOpConsolidator;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class TargetRegions implements CobaltScope
{
    private class EnrichmentMap
    {
        private final Map<HumanChromosome, ArrayList<TargetRegionEnrichment>> mChrEnrichments;

        public EnrichmentMap(
                final ListMultimap<HumanChromosome, TargetRegionEnrichment> enrichments, final Map<Chromosome,Integer> chrLengths)
        {
            mChrEnrichments = new HashMap<>();

            for(HumanChromosome chromosome : enrichments.keySet())
            {
                int length = chrLengths.get(chromosome);
                int numberOfSlots = length / WINDOW_SIZE;
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
                    position += WINDOW_SIZE;
                }

                mChrEnrichments.put(chromosome, enrichmentsForChromosome);
            }
        }

        public TargetRegionEnrichment getEnrichment(final HumanChromosome chromosome, final int position)
        {
            return mChrEnrichments.containsKey(chromosome) ? mChrEnrichments.get(chromosome).get(position / WINDOW_SIZE) : null;
        }
    }

    private final List<EnrichmentMap> mEnrichmentMaps;

    public TargetRegions()
    {
        mEnrichmentMaps = Lists.newArrayList();
    }

    public void loadNormalisationFiles(final List<String> filenames, final RefGenomeVersion refGenomeVersion)
    {
        RefGenomeCoordinates refGenomeCoordinates = RefGenomeCoordinates.refGenomeCoordinates(refGenomeVersion);

        for(String filename : filenames)
        {
            ListMultimap<HumanChromosome, TargetRegionEnrichment> chrEnrichmentMap = TargetRegionEnrichment.loadEnrichmentFile(filename);

            mEnrichmentMaps.add(new EnrichmentMap(chrEnrichmentMap, refGenomeCoordinates.Lengths));
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
    public boolean onTarget(final HumanChromosome chromosome, final int position)
    {
        // return getEnrichment(chromosome, position) != null;

        for(EnrichmentMap enrichmentMap : mEnrichmentMaps)
        {
            TargetRegionEnrichment enrichment = enrichmentMap.getEnrichment(chromosome, position);

            if(enrichment != null)
                return true;
        }

        return false;
    }

    @Override
    public double findRegionEnrichment(final HumanChromosome chromosome, final int position)
    {
        // TargetRegionEnrichment enrichment = getEnrichment(chromosome, position);
        // return enrichment == null ? -1.0 : enrichment.Enrichment;

        List<TargetRegionEnrichment> regionEnrichments = null;

        for(EnrichmentMap enrichmentMap : mEnrichmentMaps)
        {
            TargetRegionEnrichment enrichment = enrichmentMap.getEnrichment(chromosome, position);

            if(enrichment == null || !enrichment.isValid())
                continue;

            if(regionEnrichments == null)
            {
                regionEnrichments = Lists.newArrayListWithCapacity(mEnrichmentMaps.size());
            }

            regionEnrichments.add(enrichment);
        }

        if(regionEnrichments == null)
            return -1;

        // average the factors until an alternative approach is found
        return regionEnrichments.stream().mapToDouble(x -> x.Enrichment).average().orElse(0);
    }

    @VisibleForTesting
    public void addNormalisationMap(
            final ListMultimap<HumanChromosome, TargetRegionEnrichment> chrEnrichmentMap, final Map<Chromosome,Integer> chrLengths)
    {
        mEnrichmentMaps.add(new EnrichmentMap(chrEnrichmentMap, chrLengths));
    }
}
