package com.hartwig.hmftools.serve.vicc.copynumber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnownAmplificationDeletion {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract CopyNumberType type();
}
