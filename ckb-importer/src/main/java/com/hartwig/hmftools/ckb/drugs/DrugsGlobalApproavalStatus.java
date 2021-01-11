package com.hartwig.hmftools.ckb.drugs;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugsGlobalApproavalStatus {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract DrugsTherapy therapy();

    @NotNull
    public abstract DrugsIndications indications();

    @NotNull
    public abstract DrugsMolecularProfile molecularProfile();

    @NotNull
    public abstract String approvalAuthority();

    @NotNull
    public abstract String approvalStatus();

}
