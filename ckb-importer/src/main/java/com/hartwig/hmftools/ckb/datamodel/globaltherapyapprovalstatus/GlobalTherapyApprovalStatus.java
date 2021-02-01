package com.hartwig.hmftools.ckb.datamodel.globaltherapyapprovalstatus;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.GlobalApprovalStatusInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GlobalTherapyApprovalStatus {

    public abstract int totalCount();

    @NotNull
    public abstract List<GlobalApprovalStatusInfo> globalApprovalStatus();
}
