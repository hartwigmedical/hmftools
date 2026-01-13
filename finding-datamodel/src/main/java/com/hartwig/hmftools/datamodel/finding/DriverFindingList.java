package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.driver.DriverSource;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface DriverFindingList<T extends Driver> extends FindingList<T> {

    @NotNull
    default FindingsQuery<T> query() {
        return new FindingsQuery<>(all());
    }

    @NotNull
    default List<T> germlineOnly() {
        return query().driverSources(DriverSource.GERMLINE).results();
    }

    default List<T> somaticOnly() {
        return query().driverSources(DriverSource.SOMATIC).results();
    }
}
