package com.hartwig.hmftools.common.hla;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class HlaAllele {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String germlineAllele();
}