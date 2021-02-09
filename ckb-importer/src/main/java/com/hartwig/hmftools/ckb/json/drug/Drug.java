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
public abstract class Drug implements CkbJsonObject {

    public abstract int id();

    @NotNull
    public abstract String drugName();

    @NotNull
    public abstract List<String> term();

    @NotNull
    public abstract List<String> synonym();

    @Nullable
    public abstract String tradeName();

    @NotNull
    public abstract List<DescriptionInfo> description();

    @NotNull
    public abstract List<DrugClassInfo> drugClass();

    @Nullable
    public abstract String casRegistryNum();

    @Nullable
    public abstract String nctId();

    @Nullable
    public abstract Date createDate();

    @NotNull
    public abstract List<ClinicalTrialInfo> clinicalTrial();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<TherapyInfo> therapy();

    @Nullable
    public abstract List<GlobalApprovalStatusInfo> globalApprovalStatus();
}
