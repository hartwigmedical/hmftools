package com.hartwig.hmftools.ckb.drugs;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugsEvidenceReferences {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String pubMedId();

    @NotNull
    public abstract String title();

    @NotNull
    public abstract String url();
}
