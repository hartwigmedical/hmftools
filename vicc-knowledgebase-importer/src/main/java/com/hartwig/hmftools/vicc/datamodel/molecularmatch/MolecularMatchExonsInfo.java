package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchExonsInfo {

    @NotNull
    public abstract String chr();

    @NotNull
    public abstract String transcript();

    @Nullable
    public abstract String txStart();

    @Nullable
    public abstract String txEnd();

    @Nullable
    public abstract String cdsStart();

    @Nullable
    public abstract String cdsEnd();

    @NotNull
    public abstract MolecularMatchExonBoundaries exonBoundaries();
}
