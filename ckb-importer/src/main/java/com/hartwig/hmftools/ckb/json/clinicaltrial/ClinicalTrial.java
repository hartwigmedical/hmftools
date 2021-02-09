package com.hartwig.hmftools.ckb.json.clinicaltrial;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.json.CkbJsonObject;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ClinicalTrial implements CkbJsonObject {

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

    @Nullable
    public abstract Date updateDate();

    @NotNull
    public abstract List<IndicationInfo> indications();

    @NotNull
    public abstract List<ClinicalTrialVariantRequirementDetail> variantRequirementDetails();

    @NotNull
    public abstract List<ClinicalTrialLocation> clinicalTrialLocations();
}
