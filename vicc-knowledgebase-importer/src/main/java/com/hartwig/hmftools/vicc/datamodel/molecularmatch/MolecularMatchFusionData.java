package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchFusionData {

    @Nullable
    public abstract List<MolecularMatchBreg> Bgreg();

    @Nullable
    public abstract List<String> Bchr();

    @Nullable
    public abstract String synonym();

    @Nullable
    public abstract List<String> Agene();

    @Nullable
    public abstract List<String> Btx();

    @Nullable
    public abstract List<String> Achr();

    @Nullable
    public abstract List<String> ins();

    @Nullable
    public abstract String source();

    @Nullable
    public abstract List<MolecularMatchAgreg> Agreg();

    @Nullable
    public abstract List<String> Bgene();

    @Nullable
    public abstract List<String> Acoord();

    @Nullable
    public abstract List<String> Bori();

    @Nullable
    public abstract List<String> Aband();

    @Nullable
    public abstract List<String> Bband();

    @Nullable
    public abstract List<String> Aori();

    @Nullable
    public abstract List<String> Atx();

    @Nullable
    public abstract List<String> Bcoord();

    @Nullable
    public abstract String Paper();
}
