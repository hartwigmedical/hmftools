package com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsIntervention {

    @Nullable
    public abstract String interventionName();

    @NotNull
    public abstract List<String> otherNames();

    @Nullable
    public abstract String interventionType();

    @NotNull
    public abstract List<String> armGroupLabels();

    @Nullable
    public abstract String description();
}
