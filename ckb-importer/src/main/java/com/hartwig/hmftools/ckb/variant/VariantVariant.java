package com.hartwig.hmftools.ckb.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantVariant {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String fullName();

    @NotNull
    public abstract String impact();

    @NotNull
    public abstract String proteinEffect();
}
