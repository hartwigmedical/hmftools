package com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrials implements KbSpecificObject {

    @NotNull
    public abstract String status();

    @Nullable
    public abstract String startDate();

    @NotNull
    public abstract String title();

    @Nullable
    public abstract String briefTitle();

    @NotNull
    public abstract String studyType();

    @NotNull
    public abstract List<String> molecularAlterations();

    @NotNull
    public abstract String score();

    @NotNull
    public abstract List<MolecularMatchTrialsIntervention> interventions();

    @NotNull
    public abstract List<MolecularMatchTrialsLocation> locations();

    @Nullable
    public abstract MolecularMatchTrialsOverallContact overallContact();

    @NotNull
    public abstract String link();

    @NotNull
    public abstract String phase();

    @NotNull
    public abstract List<MolecularMatchTrialsTag> tags();

    @NotNull
    public abstract String id();

}
