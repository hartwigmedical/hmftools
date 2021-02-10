package com.hartwig.hmftools.ckb.datamodelinterpretation.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantInfo {

    public abstract int id();

    @NotNull
    public abstract String fullName();

    @Nullable
    public abstract String impact();

    @Nullable
    public abstract String proteinEffect();
}
