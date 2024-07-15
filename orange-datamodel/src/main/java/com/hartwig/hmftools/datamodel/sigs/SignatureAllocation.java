package com.hartwig.hmftools.datamodel.sigs;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SignatureAllocation
{
    @NotNull
    String signature();

    @NotNull
    String etiology();

    double allocation();

    double percent();
}
