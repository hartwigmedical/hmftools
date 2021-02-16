package com.hartwig.hmftools.ckb.datamodel.therapy;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GlobalTherapyApprovalStatus {

    public abstract int id();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract String approvalAuthority();

}