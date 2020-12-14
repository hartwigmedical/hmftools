package com.hartwig.hmftools.common.lims.cohort;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class LimsCohortConfigData {

    @NotNull
    public abstract String cohortId();

    @NotNull
    public abstract String reportGermline();

    @NotNull
    public abstract String reportGermlineFlag();

    @NotNull
    public abstract String reportConclusion();

    @NotNull
    public abstract String reportViral();

    @NotNull
    public abstract String requireHospitalId();

    @NotNull
    public abstract String requireHospitalPAId();
}
