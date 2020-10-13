package com.hartwig.hmftools.common.clinical;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PatientTumorLocationV2 {

    @NotNull
    public abstract String patientIdentifier();

    @NotNull
    public abstract String primaryTumorLocation();

    @NotNull
    public abstract String primaryTumorSubLocation();

    @NotNull
    public abstract String primaryTumorType();

    @NotNull
    public abstract String primaryTumorSubType();

    @NotNull
    public abstract List<String> doids();

    public abstract boolean isOverridden();

}
