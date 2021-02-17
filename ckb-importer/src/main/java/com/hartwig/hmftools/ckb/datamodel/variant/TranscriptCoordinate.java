package com.hartwig.hmftools.ckb.datamodel.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TranscriptCoordinate {

    public abstract int id();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String gDna();

    @NotNull
    public abstract String cDna();

    @NotNull
    public abstract String protein();

    @NotNull
    public abstract String sourceDb();

    @NotNull
    public abstract String refGenomeBuild();
}