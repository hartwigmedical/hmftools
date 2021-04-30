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
    public abstract String positionGRCh37();

    @NotNull
    public abstract String refGRCh37();

    @NotNull
    public abstract String altGRCh37();

    @NotNull
    public abstract String positionGRCh38();

    @NotNull
    public abstract String refGRCh38();

    @NotNull
    public abstract String altGRCh38();

    @NotNull
    public abstract String rsid();

    @NotNull
    public abstract String variantAnnotation();

    @NotNull
    public abstract String filter();

}
