package com.hartwig.hmftools.wisp.purity.loh;

import static java.lang.Math.round;

public class SiteData
{
    public final int Support;
    public final double TumorBaf;
    public final double SampleBaf;
    public final int TumorDepth;
    public final int SampleDepth;
    public final double AF;

    public SiteData(final double tumorBaf, final double sampleBaf, final int tumorDepth, final int sampleDepth)
    {
        TumorBaf = tumorBaf;
        SampleBaf = sampleBaf;
        TumorDepth = tumorDepth;
        SampleDepth = sampleDepth;

        double adjustedSampleBaf = TumorBaf < 0.5 ? SampleBaf : 1 - SampleBaf;
        double support = sampleDepth * adjustedSampleBaf;
        Support = (int) round(support);
        AF = SampleDepth > 0 ? support / SampleDepth : 0;
    }
}
