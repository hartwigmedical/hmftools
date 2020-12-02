package com.hartwig.hmftools.statcalcs.cooc;

import java.util.List;

import com.google.common.collect.Lists;

public class SampleGeneData
{
    final String SampleId;

    final List<String> GeneKnown;
    final List<String> GeneUnclear;
    final List<String> CategoryKnown;
    final List<String> CategoryUnclear;

    public SampleGeneData(final String sampleId)
    {
        SampleId = sampleId;
        GeneKnown = Lists.newArrayList();
        GeneUnclear = Lists.newArrayList();
        CategoryKnown = Lists.newArrayList();
        CategoryUnclear = Lists.newArrayList();
    }

}
