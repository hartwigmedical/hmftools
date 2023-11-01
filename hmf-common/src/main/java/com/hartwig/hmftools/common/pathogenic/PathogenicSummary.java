package com.hartwig.hmftools.common.pathogenic;

public class PathogenicSummary
{
    public final String ClinvarInfo;
    public final Pathogenicity Status;

    public PathogenicSummary(final String clinvarInfo, final Pathogenicity status)
    {
        ClinvarInfo = clinvarInfo;
        Status = status;
    }
}
