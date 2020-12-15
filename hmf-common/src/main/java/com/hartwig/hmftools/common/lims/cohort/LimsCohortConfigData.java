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

    public abstract boolean hospitalId();

    public abstract boolean reportGermline();

    public abstract boolean reportGermlineFlag();

    public abstract boolean reportConclusion();

    public abstract boolean reportViral();

    public abstract boolean requireHospitalId();

    public abstract boolean requireHospitalPAId();

    public abstract boolean hospitalPersonsStudy();

    public abstract boolean hospitalPersonsRequester();

    public abstract boolean outputFile();

    public abstract boolean submission();

    public abstract boolean sidePanelInfo();
}
