package com.hartwig.hmftools.common.hla;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LilacReporting {

    @NotNull
    public abstract LilacGermlineAllele lilacGermlineAllele();

    public abstract double germlineCopies();

    public abstract double tumorCopies();

    @NotNull
    public abstract String somaticMutations();

    @NotNull
    public abstract String interpretation();
}