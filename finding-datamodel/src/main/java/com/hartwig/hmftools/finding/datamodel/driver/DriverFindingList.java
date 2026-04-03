package com.hartwig.hmftools.finding.datamodel.driver;

import java.util.List;

import com.hartwig.hmftools.finding.datamodel.RecordBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingsQuery;
import com.hartwig.hmftools.finding.datamodel.finding.IFindingList;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record DriverFindingList<T extends Driver>(
        @NotNull FindingStatus status,
        @NotNull List<T> findings
) implements IFindingList<T>
{
    @NotNull
    public FindingsQuery<T> query()
    {
        return new FindingsQuery<>(findings());
    }

    @NotNull
    public DriverFindingList<T> germlineOnly()
    {
        return DriverFindingListBuilder.builder(this)
                .findings(query().driverSources(DriverSource.GERMLINE).results())
                .build();
    }

    @NotNull
    public DriverFindingList<T> somaticOnly()
    {
        return DriverFindingListBuilder.builder(this)
                .findings(query().driverSources(DriverSource.SOMATIC).results())
                .build();
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
        return DriverFindingListBuilder.builder(this)
                .findings(query().reportedStatuses(reportedStatuses).results())
                .build();
    }
}
