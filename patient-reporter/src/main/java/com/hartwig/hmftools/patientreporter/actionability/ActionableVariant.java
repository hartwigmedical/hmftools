package com.hartwig.hmftools.patientreporter.actionability;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableVariant {

    @NotNull
    abstract String gene();

    @NotNull
    abstract String chromosome();
}
