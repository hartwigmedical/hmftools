package com.hartwig.hmftools.cobalt.metrics;

public record WindowStatistics(String chromosome,
                               int position,
                               long count,
                               double mean,
                               double median,
                               double min,
                               double max,
                               double sd)
{
}
