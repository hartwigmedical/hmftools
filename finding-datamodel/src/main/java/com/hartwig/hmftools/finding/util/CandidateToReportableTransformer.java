package com.hartwig.hmftools.finding.util;

import static com.hartwig.hmftools.finding.util.TransformUtil.transformDriverFindingList;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.FusionBuilder;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.SmallVariantBuilder;
import com.hartwig.hmftools.finding.datamodel.Virus;
import com.hartwig.hmftools.finding.datamodel.VirusBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFields;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.ReportedStatus;

import jakarta.validation.constraints.NotNull;

public class CandidateToReportableTransformer
{
    @NotNull
    public static FindingRecord transform(@NotNull FindingRecord record)
    {
        return FindingRecordBuilder.builder(record)
                .somaticSmallVariants(transformDriverFindingList(record.somaticSmallVariants(), CandidateToReportableTransformer::transform))
                .germlineSmallVariants(transformDriverFindingList(record.germlineSmallVariants(), CandidateToReportableTransformer::transform))
                .fusions(transformDriverFindingList(record.fusions(), CandidateToReportableTransformer::transform))
                .viruses(transformDriverFindingList(record.viruses(), CandidateToReportableTransformer::transform))
                .build();
    }

    private static SmallVariant transform(SmallVariant smallVariant)
    {
        return SmallVariantBuilder.builder(smallVariant)
                .driver(transform(smallVariant.driver()))
                .build();
    }

    private static Fusion transform(Fusion fusion)
    {
        return FusionBuilder.builder(fusion)
                .driver(transform(fusion.driver()))
                .build();
    }

    private static Virus transform(Virus virus)
    {
        return VirusBuilder.builder(virus)
                .driver(transform(virus.driver()))
                .build();
    }

    private static DriverFields transform(DriverFields driverFields)
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
