package com.hartwig.hmftools.common.ecrf.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ValidationFinding {
    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String ecrfItem();

    @NotNull
    public abstract String message();

    @NotNull
    public abstract String formStatus();

    @NotNull
    public abstract String formLocked();
}
