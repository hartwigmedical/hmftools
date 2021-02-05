package com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial;

import java.util.Date;
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
    public abstract List<String> ageGroups();

    @NotNull
    public abstract String gender();

    @NotNull
    public abstract String variantRequirement();

    @NotNull
    public abstract String sponsor();

    @NotNull
    public abstract Date updateDate();

    @NotNull
    public abstract List<ClinicalTrialVariantRequirementDetail> clinicalTrialVariantRequirementDetails();

    @NotNull
    public abstract List<ClinicalTrialLocation> locations();

}



