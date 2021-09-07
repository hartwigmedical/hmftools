package com.hartwig.hmftools.neo.cohort;

public class StatusResults
{
    public double AffinityTotal;
    public int AffinityLowCount;
    public int AffinityMediumCount;

    public double PresentationTotal;
    public int PresentationCount;

    public static final int NORMAL = 0;
    public static final int TUMOR = 1;
    public static final int SIM_TUMOR = 2;
    public static final int STATUS_MAX = 3;

    public StatusResults()
    {
        AffinityLowCount = 0;
        AffinityMediumCount = 0;
        AffinityTotal = 0;
        PresentationTotal = 0;
        PresentationCount = 0;
    }

}
