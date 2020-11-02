package com.hartwig.hmftools.protect.variants.germline;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GermlineReporting {

    public abstract boolean notifyClinicalGeneticus();

    public abstract boolean reportBiallelicOnly();

    @NotNull
    public abstract String variant();
}
