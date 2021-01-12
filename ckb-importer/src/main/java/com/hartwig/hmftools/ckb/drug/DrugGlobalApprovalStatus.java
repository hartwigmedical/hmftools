package com.hartwig.hmftools.ckb.drug;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugGlobalApprovalStatus {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract DrugTherapy therapy();

    @NotNull
    public abstract DrugIndication indications();

    @NotNull
    public abstract DrugMolecularProfile molecularProfile();

    @NotNull
    public abstract String approvalAuthority();

    @NotNull
    public abstract String approvalStatus();

}
