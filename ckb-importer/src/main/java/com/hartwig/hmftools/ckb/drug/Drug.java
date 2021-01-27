package com.hartwig.hmftools.ckb.drug;

import java.util.List;

import com.hartwig.hmftools.ckb.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.common.TherapyInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Drug {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String drugName();

    @NotNull
    public abstract List<String> terms();

    @NotNull
    public abstract List<String> synonyms();

    @Nullable
    public abstract String tradeName();

    @NotNull
    public abstract List<DescriptionInfo> drugDescriptions();

    @NotNull
    public abstract List<DrugClassInfo> drugClasses();

    @Nullable
    public abstract String casRegistryNum();

    @Nullable
    public abstract String nctId();

    @NotNull
    public abstract String createDate();

    @NotNull
    public abstract List<ClinicalTrialInfo> clinicalTrials();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<TherapyInfo> therapies();

    @Nullable
    public abstract List<GlobalApprovalStatusInfo> globalApprovaStatus();


}
