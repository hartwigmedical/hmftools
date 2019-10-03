package com.hartwig.hmftools.vicc.datamodel.civic;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicClinicalTrial {

    @NotNull
    public abstract String nct_id();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract String clinical_trial_url();

    @NotNull
    public abstract String name();
}
