package com.hartwig.hmftools.finding.util;

import java.util.Objects;
import java.util.function.Function;

import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFields;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.finding.datamodel.finding.Finding;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;

import jakarta.validation.constraints.NotNull;

public class TransformUtil
{
    @NotNull
    static <T extends Driver> DriverFindingList<T> transformDriverFindingList(@NotNull DriverFindingList<T> findingList,
            Function<T, T> transformFunction)
    {
        return DriverFindingListBuilder.builder(findingList)
                .findings(findingList.stream().map(transformFunction).filter(Objects::nonNull).toList())
                .build();
    }

    @NotNull
    static <T extends Finding> FindingList<T> transformFindingList(@NotNull FindingList<T> findingList,
            Function<T, T> transformFunction)
    {
        return FindingListBuilder.builder(findingList)
                .findings(findingList.stream().map(transformFunction).filter(Objects::nonNull).toList())
                .build();
    }

    @NotNull
    static DriverFieldsBuilder toCandidate(@NotNull DriverFields driverFields)
    {
        return DriverFieldsBuilder.builder(driverFields).reportedStatus(ReportedStatus.CANDIDATE);
    }
}
