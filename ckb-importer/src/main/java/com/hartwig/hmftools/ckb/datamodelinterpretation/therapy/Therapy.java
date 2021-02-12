package com.hartwig.hmftools.ckb.datamodelinterpretation.therapy;

import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodelinterpretation.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.interpretation.common.druginterpretation.DrugInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Therapy {

    public abstract int id();

    @NotNull
    public abstract String therapyName();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract List<TherapyDescription> descriptions();

    @Nullable
    public abstract Date createDate();

    @Nullable
    public abstract Date updateDate();

    @NotNull
    public abstract DrugInterpretation drugInterpretation();

    @NotNull
    public abstract List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatuses();
}
