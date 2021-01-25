package com.hartwig.hmftools.ckb.clinicaltrial;

import java.util.List;

import com.hartwig.hmftools.ckb.common.TherapyInfo;

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
    public abstract List<TherapyInfo> therapies();

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
    public abstract List<ClinicalTrialIndication> indications();

    @NotNull
    public abstract List<ClinicalTrialVariantRequirementDetail> variantRequirementDetails();

    @NotNull
    public abstract List<ClinicalTrialLocation> clinicalTrialLocations();
}
