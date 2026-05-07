package com.hartwig.hmftools.finding.util;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;

import jakarta.validation.constraints.NotNull;

public class ReportedOnlyConverter
{
    @NotNull
    public static FindingRecord convert(@NotNull FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .somaticSmallVariants(convert(record.somaticSmallVariants()))
                .germlineSmallVariants(convert(record.germlineSmallVariants()))
                .somaticDisruptions(convert(record.somaticDisruptions()))
                .germlineDisruptions(convert(record.germlineDisruptions()))
                .somaticGainDeletions(convert(record.somaticGainDeletions()))
                .germlineGainDeletions(convert(record.germlineGainDeletions()))
                .fusions(convert(record.fusions()))
                .viruses(convert(record.viruses()))
                .build();
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> convert(@NotNull DriverFindingList<T> driverFindingList)
    {
        return driverFindingList.reportedOnly();
    }
}
