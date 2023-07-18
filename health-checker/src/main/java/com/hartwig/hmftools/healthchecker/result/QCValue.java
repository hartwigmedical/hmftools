package com.hartwig.hmftools.healthchecker.result;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class QCValue
{
    @NotNull
    public abstract QCValueType type();

    @NotNull
    public abstract String value();
}
