package com.hartwig.hmftools.common.doid;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DoidBasicPropertyValue
{
    @NotNull
    public abstract String pred();

    @NotNull
    public abstract String val();
}
