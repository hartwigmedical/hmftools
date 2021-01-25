package com.hartwig.hmftools.ckb.globaltherapyapprovalstatus;

import com.hartwig.hmftools.ckb.common.TherapyInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GlobalTherapyApprovalStatusList {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract TherapyInfo therapy();

    @NotNull
    public abstract GlobalTherapyApprovalStatusIndication indication();

    @NotNull
    public abstract GlobalTherapyApprovalStatusMolecularProfile molecularProfile();

    @NotNull
    public abstract String approvalAuthority();

    @NotNull
    public abstract String approvalStatus();
}
