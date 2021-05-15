package com.hartwig.hmftools.common.lims.cohort;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class LimsCohortConfig {

    @NotNull
    public abstract String cohortId();

    public abstract boolean sampleContainsHospitalCenterId();

    public abstract boolean reportGermline();

    public abstract boolean reportGermlineFlag();

    public abstract boolean reportConclusion();

    public abstract boolean reportViral();

    public abstract boolean reportPeach();

    public abstract boolean requireHospitalId();

    public abstract boolean requireHospitalPAId();

    public abstract boolean requireHospitalPersonsStudy();

    public abstract boolean requireHospitalPersonsRequester();

    public abstract boolean requireAdditionalInformationForSidePanel();
}
