package com.hartwig.hmftools.vicc.datamodel.jaxtrials;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JaxTrialsMolecularProfile {

    @NotNull
    public abstract String requirementType();

    @NotNull
    public abstract String profileName();

    @NotNull
    public abstract String id();
}
