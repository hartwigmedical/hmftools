package com.hartwig.hmftools.common.peach;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PeachCalls {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String positionV37();

    @NotNull
    public abstract String positionV38();

    @NotNull
    public abstract String refV37();

    @NotNull
    public abstract String refV38();

    @NotNull
    public abstract String allele1();

    @NotNull
    public abstract String allele2();

    @NotNull
    public abstract String rsid();

    @NotNull
    public abstract String variantAnnotationV37();

    @NotNull
    public abstract String filterV37();

    @NotNull
    public abstract String variantAnnotationV38();

    @NotNull
    public abstract String filterV38();

    @NotNull
    public abstract String panelVersion();

    @NotNull
    public abstract String repoVersion();
}
