package com.hartwig.hmftools.ckb.datamodel.clinicaltrial;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.indication.Indication;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrial {

    @Nullable
    public abstract Date updateDate();

    @NotNull
    public abstract String nctId();

    @NotNull
    public abstract String title();

    @NotNull
    public abstract List<Therapy> therapies();

    @NotNull
    public abstract List<Indication> indications();

    @NotNull
    public abstract String phase();

    @NotNull
    public abstract String recruitment();

    @NotNull
    public abstract List<String> ageGroups();

    @Nullable
    public abstract String gender();

    @Nullable
    public abstract String sponsor();

    @NotNull
    public abstract String variantRequirement();

    @NotNull
    public abstract List<ClinicalTrialVariantRequirementDetail> variantRequirementDetails();

    @NotNull
    public abstract List<ClinicalTrialLocation> locations();
}