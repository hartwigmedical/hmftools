package com.hartwig.hmftools.purple.fitting;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.purple.region.FittingRegion;

public class AmberPointsProvider
{
    private final GenomePositionSelector<AmberBAF> BafSelector;

    public AmberPointsProvider(final ListMultimap<Chromosome, AmberBAF> chromosomeBafs)
    {
        BafSelector = GenomePositionSelectorFactory.create(chromosomeBafs);
    }

    public List<AmberBAF> nextBlock(FittingRegion segment)
    {
        List<AmberBAF> result = new ArrayList<>();
        BafSelector.select(segment, result::add);
        return result;
    }
}
