package com.hartwig.hmftools.cobalt.calculations;

import java.util.Iterator;
import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class BamRatios
{
    ListMultimap<Chromosome, BamRatio> Ratios;

    public BamRatios(ListMultimap<Chromosome, BamRatio> ratios)
    {
        Ratios = ratios;
    }

    public void normalise(ResultsNormaliser normaliser)
    {
        Ratios.forEach(((chromosome, bamRatio) -> normaliser.recordValue(bamRatio)));
        normaliser.dataCollectionFinished();
        Ratios.forEach(((chromosome, bamRatio) -> normaliser.normalise(bamRatio)));
    }

    public void consolidate(ResultsConsolidator consolidator)
    {
        if (Ratios.isEmpty())
        {
            return;
        }
        ListMultimap<Chromosome, BamRatio> consolidatedRatios = consolidator.consolidate(Ratios);
        Preconditions.checkState(Ratios.size() >= consolidatedRatios.size());
        Preconditions.checkState(Ratios.keySet().equals(consolidatedRatios.keySet()));
        if (consolidatedRatios.size() == Ratios.size())
        {
            return;
        }

        Ratios.keySet().forEach(chromosome ->
        {
            List<BamRatio> originalRatiosForChromosome = Ratios.get(chromosome);
            List<BamRatio> consolidatedRatiosForChromosome = consolidatedRatios.get(chromosome);
            Preconditions.checkState(originalRatiosForChromosome.size() > consolidatedRatiosForChromosome.size());
            Iterator<BamRatio> originalsIterator = originalRatiosForChromosome.iterator();
            consolidatedRatiosForChromosome.forEach(consolidatedRatio ->
            {
                boolean seekingMatch = true;
                while (seekingMatch && originalsIterator.hasNext())
                {
                    BamRatio originalRatio = originalsIterator.next();
                    if (originalRatio.Position == consolidatedRatio.Position)
                    {
                        originalRatio.setRatio(consolidatedRatio.ratio());
                        seekingMatch = false;
                    }
                    else
                    {
                        originalRatio.setRatio(-1.0);
                    }
                }
            });
        });
    }

    public ListMultimap<Chromosome, BamRatio> ratios()
    {
        return ArrayListMultimap.create(Ratios);
    }
}
