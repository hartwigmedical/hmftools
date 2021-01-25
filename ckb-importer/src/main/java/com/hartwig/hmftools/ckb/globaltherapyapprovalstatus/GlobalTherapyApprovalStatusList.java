package com.hartwig.hmftools.ckb.globaltherapyapprovalstatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GlobalTherapyApprovalStatusList {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract GlobalTherapyApprovalStatusTherapy therapy();

    @NotNull
    public abstract GlobalTherapyApprovalStatusIndication indication();

    @NotNull
    public abstract GlobalTherapyApprovalStatusMolecularProfile molecularProfile();

    @NotNull
    public abstract String approvalAuthority();

    @NotNull
    public abstract String approvalStatus();
}
