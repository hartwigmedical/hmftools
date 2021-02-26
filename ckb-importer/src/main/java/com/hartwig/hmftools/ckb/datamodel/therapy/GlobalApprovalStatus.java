package com.hartwig.hmftools.ckb.datamodel.therapy;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GlobalApprovalStatus {

    public abstract int id();

    public abstract int profileId();

    public abstract int indicationId();

    @NotNull
    public abstract String approvalStatus();

    @NotNull
    public abstract String approvalAuthority();

}