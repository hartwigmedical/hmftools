package com.hartwig.hmftools.common.actionability.somaticvariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableSomaticVariant {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String chromosome();

    public abstract long position();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract String reference();

    @NotNull
    public abstract String drug();

    @NotNull
    public abstract String drugsType();

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract String level();

    @NotNull
    public abstract String response();
}
