package com.hartwig.hmftools.finding;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.finding.ReportedStatus;

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
}
