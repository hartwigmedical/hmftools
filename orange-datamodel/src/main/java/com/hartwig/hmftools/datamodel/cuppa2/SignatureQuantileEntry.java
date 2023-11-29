package com.hartwig.hmftools.datamodel.cuppa2;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SignatureQuantileEntry
{
    String featName();
    double featValue();
    String cancerType();
    double dataValue();
    int rank();
    int rankGroup();
}

