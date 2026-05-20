package com.hartwig.hmftools.finding.util;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;

import jakarta.validation.constraints.NotNull;

public class ReportedOnlyTransformer
{
    @NotNull
    public static FindingRecord transform(@NotNull FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .somaticSmallVariants(transform(record.somaticSmallVariants()))
                .germlineSmallVariants(transform(record.germlineSmallVariants()))
                .somaticDisruptions(transform(record.somaticDisruptions()))
                .germlineDisruptions(transform(record.germlineDisruptions()))
                .somaticGainDeletions(record.somaticGainDeletions())
                .germlineGainDeletions(transform(record.germlineGainDeletions()))
                .fusions(transform(record.fusions()))
                .viruses(transform(record.viruses()))
                .build();
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> transform(@NotNull DriverFindingList<T> driverFindingList)
    {
        return driverFindingList.reportedOnly();
    }
}
