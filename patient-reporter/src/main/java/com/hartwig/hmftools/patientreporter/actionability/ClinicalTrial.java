package com.hartwig.hmftools.patientreporter.actionability;

import com.hartwig.hmftools.common.actionability.ActionabilitySource;

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
    public abstract String level();

    @NotNull
    public abstract String acronym();

    @NotNull
    public abstract ActionabilitySource source();

    @NotNull
    public abstract String reference();

    public abstract boolean isOnLabel();
}
