package com.hartwig.hmftools.common.actionability.cancertype;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class CancerTypeReading {

    @NotNull
    abstract String cancerType();

    @NotNull
    abstract String doidSet();
}
