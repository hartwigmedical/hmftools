package com.hartwig.hmftools.viridian.reporting;

// Why a virus was or was not reported.
public enum VirusReportStatus
{
    // Reported because the virus genome is integrated into the tumor genome.
    REPORTED_ON_INTEGRATION(true),
    // Reported because the virus genome is highly present in each tumor cell.
    REPORTED_ON_CLONALITY(true),
    REPORTED_ON_INTEGRATION_AND_CLONALITY(true),
    // Not a reportable virus.
    NOT_REPORTABLE(false),
    // Reportable but no integrations and clonality couldn't be evaluated due to failing Purple data.
    CLONALITY_UNEVALUABLE(false),
    // Reportable but not enough virus genome present per tumor cell.
    NOT_CLONAL(false);

    private final boolean mIsReported;

    VirusReportStatus(boolean isReported)
    {
        mIsReported = isReported;
    }

    public boolean isReported()
    {
        return mIsReported;
    }
}
