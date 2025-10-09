package com.hartwig.hmftools.cobalt.calculations;

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

    }

    public ListMultimap<Chromosome, BamRatio> ratios()
    {
        return ArrayListMultimap.create(Ratios);
    }
}
