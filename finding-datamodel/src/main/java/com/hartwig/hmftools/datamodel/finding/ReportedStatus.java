package com.hartwig.hmftools.datamodel.finding;

public enum ReportedStatus
{
    // must be listed this way to allow comparison
    NON_DRIVER_GENE,
    NOT_REPORTED,
    CANDIDATE,
    REPORTED;

    public static boolean isMoreReportable(ReportedStatus r1, ReportedStatus r2)
    {
        return r1.ordinal() > r2.ordinal();
    }
}
