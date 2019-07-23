package com.hartwig.hmftools.linx.stats;

import java.util.List;

import com.google.common.collect.Lists;

public class SampleStatsData
{
    final String SampleId;

    final List<String> GeneKnown;
    final List<String> GeneUnclear;
    final List<String> CategoryKnown;
    final List<String> CategoryUnclear;

    public SampleStatsData(final String sampleId)
    {
        SampleId = sampleId;
        GeneKnown = Lists.newArrayList();
        GeneUnclear = Lists.newArrayList();
        CategoryKnown = Lists.newArrayList();
        CategoryUnclear = Lists.newArrayList();
    }

}
