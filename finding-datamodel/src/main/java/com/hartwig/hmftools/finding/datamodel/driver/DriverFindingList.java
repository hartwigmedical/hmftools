package com.hartwig.hmftools.finding.datamodel.driver;

import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.finding.datamodel.RecordBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.IFindingList;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record DriverFindingList<T extends Driver>(
        @NotNull FindingStatus status,
        @NotNull List<T> findings,
        @Nullable Double purityThreshold
) implements IFindingList<T>
{
    @NotNull
    public DriverFindingList<T> germlineOnly()
    {
        return filter(f -> f.driverSource() == DriverSource.GERMLINE);
    }

    @NotNull
    public DriverFindingList<T> somaticOnly()
    {
        return filter(f -> f.driverSource() == DriverSource.SOMATIC);
    }

    @NotNull
    public DriverFindingList<T> reportedOnly()
    {
        return filter(f -> f.reportedStatus() == ReportedStatus.REPORTED);
    }

    @NotNull
    public DriverFindingList<T> candidateOnly()
    {
        return filter(f -> f.reportedStatus() == ReportedStatus.CANDIDATE);
    }

    @NotNull
    public DriverFindingList<T> filter(Predicate<T> filter)
    {
        return DriverFindingListBuilder.builder(this)
                .findings(findings.stream().filter(filter).toList())
                .build();
    }

    public List<T> findingsIfOk()
    {
        return status().isOK() ? findings() : List.of();
    }
}
