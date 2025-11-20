package com.hartwig.hmftools.orange.algo.util;

import java.util.Arrays;
import java.util.Collection;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;

import org.jetbrains.annotations.NotNull;

public class DriverUtils
{
    public static ReportedStatus convertReportedStatus(com.hartwig.hmftools.common.purple.ReportedStatus reportedStatus)
    {
        return switch (reportedStatus)
        {
            case NONE -> ReportedStatus.NONE;
            case CANDIDATE -> ReportedStatus.CANDIDATE;
            case REPORTED -> ReportedStatus.REPORTED;
        };
    }

    public static ReportedStatus maxReportedStatus(@NotNull ReportedStatus... reportedStatuses)
    {
        return maxReportedStatus(Arrays.asList(reportedStatuses));
    }

    public static ReportedStatus maxReportedStatus(@NotNull Collection<ReportedStatus> reportedStatuses)
    {
        ReportedStatus maxStatus = ReportedStatus.NONE;

        for(ReportedStatus reportedStatus : reportedStatuses)
        {
            if(reportedStatus == ReportedStatus.REPORTED)
            {
                return ReportedStatus.REPORTED;
            }
            else if(reportedStatus == ReportedStatus.CANDIDATE)
            {
                maxStatus = ReportedStatus.CANDIDATE;
            }
        }
        return maxStatus;
    }

    public static DriverInterpretation maxDriverInterpretation(@NotNull DriverInterpretation... driverInterpretations)
    {
        return maxDriverInterpretation(Arrays.asList(driverInterpretations));
    }

    public static DriverInterpretation maxDriverInterpretation(@NotNull Collection<DriverInterpretation> driverInterpretations)
    {
        DriverInterpretation maxDriverInterpretation = DriverInterpretation.LOW;

        for(DriverInterpretation driverInterpretation : driverInterpretations)
        {
            if(driverInterpretation == DriverInterpretation.HIGH)
            {
                return DriverInterpretation.HIGH;
            }
            if(driverInterpretation == DriverInterpretation.MEDIUM)
            {
                maxDriverInterpretation = DriverInterpretation.MEDIUM;
            }
        }
        return maxDriverInterpretation;
    }
}
