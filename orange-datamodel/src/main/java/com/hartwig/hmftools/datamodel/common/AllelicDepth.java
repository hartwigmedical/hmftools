package com.hartwig.hmftools.datamodel.common;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface AllelicDepth
{
    int totalReadCount();
    int alleleReadCount();

    default double alleleFrequency()
    {
        return totalReadCount() > 0 ? alleleReadCount() / (double)totalReadCount() : 0;
    }
}
