package com.hartwig.hmftools.common.clinical;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PatientTumorLocation {

    @NotNull
    public abstract String patientIdentifier();

    @NotNull
    public abstract String primaryTumorLocation();

    @NotNull
    public abstract String cancerSubtype();
}
