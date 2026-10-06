package com.hartwig.hmftools.viridian.reporting;

// Why a virus was or was not reported.
public enum VirusReportStatus
{
    REPORTED_ON_INTEGRATION(true),
    REPORTED_ON_COPY_NUMBER(true),
    // Not a reportable virus.
    NOT_REPORTABLE(false),
    // Reportable but copy number could not be evaluated: an unresolved group, or an unusable Purple fit.
    COPY_NUMBER_UNEVALUABLE(false),
    // Reportable but not enough virus genome present per tumor cell.
    COPY_NUMBER_TOO_LOW(false);

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
