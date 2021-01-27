package com.hartwig.hmftools.ckb.datamodel.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantTranscriptCoordinate {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String gDNA();

    @NotNull
    public abstract String cDNA();

    @NotNull
    public abstract String protein();

    @NotNull
    public abstract String sourceDB();

    @NotNull
    public abstract String refGenomeBuild();
}
