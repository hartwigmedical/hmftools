package com.hartwig.hmftools.ckb.datamodelinterpretation.globaltherapyapprovalstatus;

import com.hartwig.hmftools.ckb.datamodelinterpretation.indication.Indication;
import com.hartwig.hmftools.ckb.interpretation.common.variantinterpretation.VariantInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GlobalTherapyApprovalStatus {

    public abstract int id();

    @NotNull
    public abstract VariantInterpretation variantInterpretation();

    @NotNull
    public abstract Indication indications();

    @NotNull
    public abstract String approvalAuthority();

    @NotNull
    public abstract String approvalStatus();
}