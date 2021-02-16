package com.hartwig.hmftools.ckb.datamodel.therapy;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GlobalTherapyApprovalStatus {

    public abstract int id();

    public abstract int profileId();

    // TODO Check if always equal to the therapy the approval status belongs to?
    public abstract int therapyId();

    @NotNull
    public abstract String indicationId();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract String approvalAuthority();

}