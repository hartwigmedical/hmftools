package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchWGSaMap {

    @Nullable
    public abstract String AA();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String GRCh37_Chr_Start_Ref_Alt();

    @NotNull
    public abstract List<String> Synonyms();

    @NotNull
    public abstract List<String> ProtCoords();

    @NotNull
    public abstract String NucleotideChange();

    @Nullable
    public abstract String Exon();

    @NotNull
    public abstract String Gene();

    @NotNull
    public abstract String Transcript();
}
