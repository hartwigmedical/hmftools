package com.hartwig.hmftools.serve.vicc;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableEvidence {

    @NotNull
    public abstract String drug();

    @NotNull
    public abstract String drugType();

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract String cancerTypeDOID();

    @NotNull
    public abstract String evidenceLevel();

    @NotNull
    public abstract String evidenceDirection();
}
