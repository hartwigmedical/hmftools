package com.hartwig.hmftools.cobalt.targeted;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.WINDOW_SIZE;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.calculations.BamRatio;
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
    private final Map<HumanChromosome, ArrayList<TargetRegionEnrichment>> mChrEnrichments;

    public TargetRegions()
    {
        mChrEnrichments = new HashMap<>();
    }

    public void loadNormalisationFile(final String filename, final RefGenomeVersion refGenomeVersion)
    {
        RefGenomeCoordinates refGenomeCoordinates = RefGenomeCoordinates.refGenomeCoordinates(refGenomeVersion);
        ListMultimap<HumanChromosome, TargetRegionEnrichment> chrEnrichmentMap = TargetRegionEnrichment.loadEnrichmentFile(filename);
        addEnrichmentData(chrEnrichmentMap, refGenomeCoordinates.Lengths);
    }

    private void addEnrichmentData(
            final ListMultimap<HumanChromosome, TargetRegionEnrichment> chrEnrichmentMap, final Map<Chromosome,Integer> chrLengths)
    {
        for(HumanChromosome chromosome : chrEnrichmentMap.keySet())
        {
            int length = chrLengths.get(chromosome);
            int numberOfSlots = length / WINDOW_SIZE;
            ArrayList<TargetRegionEnrichment> enrichmentsForChromosome = new ArrayList<>(numberOfSlots);
            int position = 1;

            Map<Integer, TargetRegionEnrichment> positionToSuppliedItem = new HashMap<>();

            for(TargetRegionEnrichment enrichment : chrEnrichmentMap.get(chromosome))
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
        return getEnrichment(chromosome, position) != null;
    }

    @Override
    public double findRegionEnrichment(final HumanChromosome chromosome, final int position)
    {
        TargetRegionEnrichment enrichment = getEnrichment(chromosome, position);
        return enrichment == null ? -1.0 : enrichment.Enrichment;
    }

    private TargetRegionEnrichment getEnrichment(final HumanChromosome chromosome, final int position)
    {
        return mChrEnrichments.containsKey(chromosome) ? mChrEnrichments.get(chromosome).get(position / WINDOW_SIZE) : null;
    }

    public static void mergeTumorRatios(
            final ListMultimap<HumanChromosome, BamRatio> allRatios, final ListMultimap<HumanChromosome, BamRatio> panelRatios)
    {
        if(allRatios.isEmpty())
        {
            allRatios.putAll(panelRatios);
            return;
        }

        // other merge ratios, removing duplicates
        for(HumanChromosome chromosome : panelRatios.keySet())
        {
            List<BamRatio> newRatios = panelRatios.get(chromosome);
            List<BamRatio> existingRatios = allRatios.get(chromosome);

            if(existingRatios == null)
            {
                allRatios.putAll(chromosome, newRatios);
                continue;
            }

            if(newRatios.size() != existingRatios.size())
            {
                CB_LOGGER.error("chromosome({}) inconsistent tumor ratio array size(existing={} new={})",
                        chromosome, existingRatios.size(), newRatios.size());
                System.exit(1);
            }

            // NOTE: if both ratios are valid then no attempt is made to reconcile the values
            for(int i = 0; i < existingRatios.size(); ++i)
            {
                BamRatio newRatio = newRatios.get(i);
                BamRatio existingRatio = existingRatios.get(i);

                if(!existingRatio.isValid() && newRatio.isValid())
                {
                    existingRatio.overrideRatio(newRatio.ratio());
                }
            }
        }
    }

    @VisibleForTesting
    public void addNormalisationMap(
            final ListMultimap<HumanChromosome, TargetRegionEnrichment> chrEnrichmentMap, final Map<Chromosome,Integer> chrLengths)
    {
        addEnrichmentData(chrEnrichmentMap, chrLengths);
    }
}
