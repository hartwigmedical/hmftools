package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrials implements KbSpecificObject {

    @NotNull
    public abstract String status();

    @NotNull
    public abstract String title();

    @NotNull
    public abstract List<String> molecularAlterations();

    @NotNull
    public abstract String score();

    @NotNull
    public abstract List<MolecularMatchTrialsIntervation> intervation();

    @NotNull
    public abstract List<MolecularMatchTrialsLocations> locations();

    @Nullable
    public abstract String briefTitle();

    @Nullable
    public abstract MolecularMatchTrialsOverallContact overallContact();

    @NotNull
    public abstract String link();

    @NotNull
    public abstract String phase();

    @NotNull
    public abstract List<MolecularMatchTrialsTags> tags();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String studyType();


}
