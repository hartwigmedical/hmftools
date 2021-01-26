package com.hartwig.hmftools.ckb.drug;

import com.hartwig.hmftools.ckb.common.IndicationInfo;
import com.hartwig.hmftools.ckb.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.common.TherapyInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugGlobalApprovalStatus {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract TherapyInfo therapy();

    @NotNull
    public abstract IndicationInfo indications();

    @NotNull
    public abstract MolecularProfileInfo molecularProfile();

    @NotNull
    public abstract String approvalAuthority();

    @NotNull
    public abstract String approvalStatus();

}
