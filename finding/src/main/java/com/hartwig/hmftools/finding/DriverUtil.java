package com.hartwig.hmftools.finding;

import com.hartwig.hmftools.finding.datamodel.DriverInterpretation;
import com.hartwig.hmftools.finding.datamodel.ReportedStatus;

@SuppressWarnings("unused")
final class DriverUtil
{
    static ReportedStatus reportedStatus(boolean isDriverGene, boolean isReportable, DriverInterpretation driverInterpretation)
    {
        if(!isDriverGene)
        {
            return ReportedStatus.NON_DRIVER_GENE;
        }
        if(isReportable)
        {
            if(driverInterpretation == DriverInterpretation.LOW)
                return ReportedStatus.CANDIDATE;
            else
                return ReportedStatus.REPORTED;
        }
        else
        {
            return ReportedStatus.NOT_REPORTED;
        }
    }

    static ReportedStatus reportedStatus(com.hartwig.hmftools.datamodel.driver.ReportedStatus purpleReportedStatus,
            DriverInterpretation driverInterpretation)
    {
        switch (purpleReportedStatus)
        {
            case NONE: {
                return ReportedStatus.NON_DRIVER_GENE;
            }
            case NOT_REPORTED: {
                return ReportedStatus.NOT_REPORTED;
            }
            case REPORTED: {
                if(driverInterpretation == DriverInterpretation.LOW)
                    return ReportedStatus.CANDIDATE;
                else
                    return ReportedStatus.REPORTED;
            }
            default: throw new IllegalArgumentException("Unhandled purple reported status: " + purpleReportedStatus);
        }
    }

    static ReportedStatus reportedStatus(com.hartwig.hmftools.datamodel.driver.ReportedStatus purpleReportedStatus)
    {
        return switch (purpleReportedStatus) {
            case NONE -> ReportedStatus.NON_DRIVER_GENE;
            case NOT_REPORTED -> ReportedStatus.NOT_REPORTED;
            case REPORTED -> ReportedStatus.REPORTED;
        };
    }

    static DriverInterpretation convert(com.hartwig.hmftools.datamodel.driver.DriverInterpretation driverInterpretation)
    {
        return switch (driverInterpretation)
        {
            case UNKNOWN -> DriverInterpretation.UNKNOWN;
            case LOW -> DriverInterpretation.LOW;
            case MEDIUM -> DriverInterpretation.MEDIUM;
            case HIGH -> DriverInterpretation.HIGH;
        };
    }
}
