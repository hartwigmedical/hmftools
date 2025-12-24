package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Findings<T extends Finding> {

    @NotNull
    FindingsStatus status();

    @NotNull
    List<T> all();
}
