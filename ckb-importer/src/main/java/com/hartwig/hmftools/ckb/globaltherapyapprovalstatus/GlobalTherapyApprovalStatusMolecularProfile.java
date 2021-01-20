package com.hartwig.hmftools.ckb.globaltherapyapprovalstatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GlobalTherapyApprovalStatusMolecularProfile {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String profileName();
}
