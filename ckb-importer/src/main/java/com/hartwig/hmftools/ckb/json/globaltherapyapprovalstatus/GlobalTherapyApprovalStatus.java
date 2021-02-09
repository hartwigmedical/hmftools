package com.hartwig.hmftools.ckb.json.globaltherapyapprovalstatus;

import java.util.List;

import com.hartwig.hmftools.ckb.json.CkbJsonObject;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GlobalTherapyApprovalStatus implements CkbJsonObject {

    public abstract int totalCount();

    @NotNull
    public abstract List<GlobalApprovalStatusInfo> globalApprovalStatus();
}
