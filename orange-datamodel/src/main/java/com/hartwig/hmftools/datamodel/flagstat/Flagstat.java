package com.hartwig.hmftools.datamodel.flagstat;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Flagstat
{
    long uniqueReadCount();

    long secondaryCount();

    long supplementaryCount();

    double mappedProportion();
}
