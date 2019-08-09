package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchExonInfo {

    @NotNull
    public abstract MolecularMatchExonsBoundries exonBoundaries();

    @Nullable
    public abstract String txStart();

    @Nullable
    public abstract String cdsEnd();

    @NotNull
    public abstract String chr();

    @Nullable
    public abstract String cdsStart();

    @NotNull
    public abstract String transcript();

    @Nullable
    public abstract String txEnd();
}
