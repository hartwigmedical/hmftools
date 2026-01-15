package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.driver.DriverSource;

import org.jetbrains.annotations.NotNull;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record DriverFindingList<T extends Driver>(
        @NotNull FindingsStatus status,
        @NotNull List<T> all
) {
    @NotNull
    public FindingsQuery<T> query() {
        return new FindingsQuery<>(all());
    }

    @NotNull
    public List<T> germlineOnly() {
        return query().driverSources(DriverSource.GERMLINE).results();
    }

    public List<T> somaticOnly() {
        return query().driverSources(DriverSource.SOMATIC).results();
    }
}
