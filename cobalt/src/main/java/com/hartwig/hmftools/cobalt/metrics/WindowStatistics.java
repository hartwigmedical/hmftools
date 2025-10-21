package com.hartwig.hmftools.cobalt.metrics;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public record WindowStatistics(String chromosome,
                               int position,
                               long count,
                               double mean,
                               double median,
                               double max,
                               double min,
                               double sd)
{
}
