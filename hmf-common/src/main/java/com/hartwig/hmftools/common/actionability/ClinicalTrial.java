package com.hartwig.hmftools.common.actionability;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrial {

    @NotNull
    public abstract String event();

    @NotNull
    public abstract String acronym();

    @NotNull
    public abstract ActionabilitySource source();

    @NotNull
    public abstract String reference();

    public abstract boolean isOnLabel();

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract EvidenceScope scope();
}
