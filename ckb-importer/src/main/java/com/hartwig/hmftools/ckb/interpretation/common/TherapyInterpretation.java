package com.hartwig.hmftools.ckb.interpretation.common;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodelinterpretation.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.datamodelinterpretation.therapy.Therapy;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TherapyInterpretation {

    @NotNull
    public abstract Therapy therapy();

    @NotNull
    public abstract List<DrugsInterpretation> drugs();

    @NotNull
    public abstract List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatuses();

}
