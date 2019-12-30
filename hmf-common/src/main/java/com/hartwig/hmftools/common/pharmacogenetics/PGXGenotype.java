package com.hartwig.hmftools.common.pharmacogenetics;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PGXGenotype {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String haplotype();

    @NotNull
    public abstract String function();

    @NotNull
    public abstract String linkedDrugs();

    @NotNull
    public abstract String urlPrescriptionInfo();

    @NotNull
    public abstract String panelVersion();

    @NotNull
    public abstract String repoVersion();
}
