package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Jax implements KbSpecificObject {

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract JaxMolecularProfile molecularProfile();

    @NotNull
    public abstract JaxTherapy therapy();

    @NotNull
    public abstract String evidenceType();

    @NotNull
    public abstract JaxIndications indications();

    @NotNull
    public abstract String efficacyEvidence();

    @Nullable
    public abstract List<JaxReferences> references();

    @NotNull
    public abstract String id();

}
