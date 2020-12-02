package com.hartwig.hmftools.vicc.datamodel.civic;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicClinicalTrial {

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String nctId();

    @NotNull
    public abstract String clinicalTrialUrl();

    @NotNull
    public abstract String description();
}
