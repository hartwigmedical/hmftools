package com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrialLocation {

    @NotNull
    public abstract String nctId();

    @NotNull
    public abstract String facility();

    @NotNull
    public abstract String city();

    @NotNull
    public abstract String country();

    @NotNull
    public abstract String status();

    @NotNull
    public abstract String state();

    @NotNull
    public abstract String zip();

    @NotNull
    public abstract List<ClinicalTrialContact> contacts();
}

