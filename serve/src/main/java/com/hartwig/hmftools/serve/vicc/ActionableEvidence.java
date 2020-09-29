package com.hartwig.hmftools.serve.vicc;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableEvidence {

    // TODO Should possibly turn into List<String> drug() and parse.
    @NotNull
    public abstract String drugs();

    @NotNull
    public abstract String cancerTypeString();

    @NotNull
    public abstract String cancerTypeDOID();

    @NotNull
    public abstract String level();

    @NotNull
    public abstract String direction();
}
