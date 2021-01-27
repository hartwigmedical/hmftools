package com.hartwig.hmftools.ckb.datamodel.drug;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.datamodel.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.datamodel.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.datamodel.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;

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

    @NotNull
    public abstract String createDate();

    @NotNull
    public abstract List<ClinicalTrialInfo> clinicalTrial();

    @NotNull
    public abstract List<EvidenceInfo> evidence();

    @NotNull
    public abstract List<TherapyInfo> therapy();

    @Nullable
    public abstract List<GlobalApprovalStatusInfo> globalApprovalStatus();


}
