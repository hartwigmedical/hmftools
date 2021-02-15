package com.hartwig.hmftools.ckb.json.drug;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.json.CkbJsonObject;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JsonDrug implements CkbJsonObject {

    public abstract int id();

    @NotNull
    public abstract String drugName();

    @NotNull
    public abstract List<String> terms();

    @NotNull
    public abstract List<String> synonyms();

    @Nullable
    public abstract String tradeName();

    @NotNull
    public abstract List<DescriptionInfo> descriptions();

    @NotNull
    public abstract List<DrugClassInfo> drugClasses();

    @Nullable
    public abstract String casRegistryNum();

    @Nullable
    public abstract String ncitId();

    @Nullable
    public abstract Date createDate();

    @NotNull
    public abstract List<ClinicalTrialInfo> clinicalTrials();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<TherapyInfo> therapies();

    @Nullable
    public abstract List<GlobalApprovalStatusInfo> globalApprovalStatus();
}
