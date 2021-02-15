package com.hartwig.hmftools.ckb.json.therapy;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.json.CkbJsonObject;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.DrugInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class JsonTherapy implements CkbJsonObject {

    public abstract int id();

    @NotNull
    public abstract String therapyName();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract List<DescriptionInfo> descriptions();

    @Nullable
    public abstract Date createDate();

    @Nullable
    public abstract Date updateDate();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<ClinicalTrialInfo> clinicalTrials();

    @NotNull
    public abstract List<DrugInfo> drugs();

    @NotNull
    public abstract List<GlobalApprovalStatusInfo> globalApprovalStatuses();
}
