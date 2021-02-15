package com.hartwig.hmftools.ckb.datamodel.therapy;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.Indication;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

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

    @NotNull
    public abstract List<Variant> variants();

    @NotNull
    public abstract Indication indication();
}