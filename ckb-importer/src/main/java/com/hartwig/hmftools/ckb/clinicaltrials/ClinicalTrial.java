package com.hartwig.hmftools.ckb.clinicaltrials;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrial {

    @NotNull
    public abstract String nctId();

    @NotNull
    public abstract String title();

    @NotNull
    public abstract String phase();

    @NotNull
    public abstract String recruitment();

    @NotNull
    public abstract List<Therapies> therapies();

    @NotNull
    public abstract List<String> ageGroups();

    @Nullable
    public abstract String gender();

    @NotNull
    public abstract String variantRequirements();

    @Nullable
    public abstract String sponsors();

    @NotNull
    public abstract String updateDate();

    @NotNull
    public abstract List<Indications> indications();

    @NotNull
    public abstract List<ClinicalTrialVariantRequirementDetails> variantRequirementDetails();

    @NotNull
    public abstract List<ClinicalTrialLocations> clinicalTrialLocations();
}
