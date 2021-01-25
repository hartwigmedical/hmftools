package com.hartwig.hmftools.ckb.globaltherapyapprovalstatus;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GlobalTherapyApprovalStatus {

    @NotNull
    public abstract String totalCount();

    @NotNull
    public abstract List<GlobalTherapyApprovalStatusList> approvalStatusList();
}
