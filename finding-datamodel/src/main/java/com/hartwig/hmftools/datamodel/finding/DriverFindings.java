package com.hartwig.hmftools.datamodel.finding;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface DriverFindings<T extends Driver> extends Findings<T> {

    @NotNull
    default List<T> allFindings(ReportedStatus... reportedStatuses)
    {
        return filtered(findings(), Set.of(), Set.of(reportedStatuses));
    }

    @NotNull
    default List<T> somaticFindings(ReportedStatus... reportedStatuses)
    {
        return filtered(findings(), Set.of(DriverSource.SOMATIC), Set.of(reportedStatuses));
    }

    @NotNull
    default List<T> germlineFindings(ReportedStatus... reportedStatuses)
    {
        return filtered(findings(), Set.of(DriverSource.GERMLINE), Set.of(reportedStatuses));
    }

    @NotNull
    private static <T extends Driver> List<T> filtered(@NotNull List<T> drivers, @NotNull Set<DriverSource> driverSources, @NotNull Set<ReportedStatus> reportedStatuses)
    {
        return drivers.stream()
                .filter(o -> driverSources.isEmpty() || driverSources.contains(o.driverSource()))
                .filter(o -> reportedStatuses.isEmpty() || reportedStatuses.contains(o.reportedStatus())).toList();
    }
}
