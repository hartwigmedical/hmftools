package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.driver.DriverSource;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record DriverFindingList<T extends Driver>(
        @NotNull FindingsStatus status,
        @NotNull List<T> all
) implements IFindingList<T>
{
    @NotNull
    public FindingsQuery<T> query()
    {
        return new FindingsQuery<>(all());
    }

    @NotNull
    public DriverFindingList<T> germlineOnly()
    {
        return new DriverFindingList<>(status, query().driverSources(DriverSource.GERMLINE).results());
    }

    @NotNull
    public DriverFindingList<T> somaticOnly()
    {
        return new DriverFindingList<>(status, query().driverSources(DriverSource.SOMATIC).results());
    }

    @NotNull
    public DriverFindingList<T> reportedOnly()
    {
        return filter(ReportedStatus.REPORTED);
    }

    @NotNull
    public DriverFindingList<T> candidateOnly()
    {
        return filter(ReportedStatus.CANDIDATE);
    }

    @NotNull
    public DriverFindingList<T> filter(ReportedStatus... reportedStatuses)
    {
        return new DriverFindingList<>(status, query().reportedStatuses(reportedStatuses).results());
    }
}
