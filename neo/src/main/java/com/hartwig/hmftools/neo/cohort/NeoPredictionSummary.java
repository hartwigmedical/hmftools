package com.hartwig.hmftools.neo.cohort;

import java.util.Set;

import com.google.common.collect.Sets;

public class NeoPredictionSummary
{
    public final Set<String> NineMers;
    public double MinAffinity;
    public double AffinityTotal;
    public double MaxPresentation;
    public double PresentationTotal;

    public NeoPredictionSummary()
    {
        NineMers = Sets.newHashSet();
        MinAffinity = 0;
        AffinityTotal = 0;
        MaxPresentation = 0;
        PresentationTotal = 0;
    }

}
