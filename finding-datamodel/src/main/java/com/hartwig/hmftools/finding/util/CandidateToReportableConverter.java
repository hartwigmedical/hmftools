package com.hartwig.hmftools.finding.util;

import java.util.function.Function;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.FusionBuilder;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.SmallVariantBuilder;
import com.hartwig.hmftools.finding.datamodel.Virus;
import com.hartwig.hmftools.finding.datamodel.VirusBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFields;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.ReportedStatus;

import jakarta.validation.constraints.NotNull;

public class CandidateToReportableConverter
{
    @NotNull
    public static FindingRecord convert(@NotNull FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .somaticSmallVariants(convert(record.somaticSmallVariants(), CandidateToReportableConverter::convert))
                .germlineSmallVariants(convert(record.germlineSmallVariants(), CandidateToReportableConverter::convert))
                .fusions(convert(record.fusions(), CandidateToReportableConverter::convert))
                .viruses(convert(record.viruses(), CandidateToReportableConverter::convert))
                .build();
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> convert(@NotNull DriverFindingList<T> driverFindingList,
            Function<T, T> convertFunction)
    {
        return DriverFindingListBuilder.builder(driverFindingList)
                .findings(driverFindingList.stream().map(convertFunction).toList())
                .build();
    }

    private static SmallVariant convert(SmallVariant smallVariant)
    {
        return SmallVariantBuilder.builder(smallVariant)
                .driver(convert(smallVariant.driver()))
                .build();
    }

    private static Fusion convert(Fusion fusion)
    {
        return FusionBuilder.builder(fusion)
                .driver(convert(fusion.driver()))
                .build();
    }

    private static Virus convert(Virus virus)
    {
        return VirusBuilder.builder(virus)
                .driver(convert(virus.driver()))
                .build();
    }

    private static DriverFields convert(DriverFields driverFields)
    {
        if(driverFields.reportedStatus() == ReportedStatus.CANDIDATE)
        {
            return DriverFieldsBuilder.builder(driverFields)
                    .reportedStatus(ReportedStatus.REPORTED)
                    .build();
        }
        return driverFields;
    }
}
