package com.hartwig.hmftools.common.purple;

public enum ReportedStatus
{
    NONE, // not used in the driver catalog, but a indicator in other types of files that a gene is not in the driver gene panel
    NOT_REPORTED, // the driver type is set to false for its reportability type in the driver gene panel
    REPORTED; // the driver type is set to true for its reportability
}
