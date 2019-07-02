package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JaxTrialsVariantRequirementDetails {

    @NotNull
    public abstract List<JaxTrialsMolecularProfile> molecularProfiles();

    @NotNull
    public abstract String requirementType();
}
