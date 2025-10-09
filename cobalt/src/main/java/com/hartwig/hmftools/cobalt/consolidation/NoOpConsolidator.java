package com.hartwig.hmftools.cobalt.consolidation;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.consolidation.ResultsConsolidator.consolidateIntoBuckets;

import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.calculations.BamRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class NoOpConsolidator implements ResultsConsolidator
{
}
