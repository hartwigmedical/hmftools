package com.hartwig.hmftools.ckb.drugs;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugsTherapies {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String therapyName();

    @NotNull
    public abstract String synonyms();
}
