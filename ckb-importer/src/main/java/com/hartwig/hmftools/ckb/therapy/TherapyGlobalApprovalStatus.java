package com.hartwig.hmftools.ckb.therapy;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TherapyGlobalApprovalStatus {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract TherapyTherapy therapy();

    @NotNull
    public abstract TherapyIndication indication();

    @NotNull
    public abstract TherapyMolecularProfile molecularProfile();

    @NotNull
    public abstract String approvalAuthority();

    @NotNull
    public abstract String approvalStatus();
}
