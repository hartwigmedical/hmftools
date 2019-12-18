package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchWGSAMap {

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcript();

    @Nullable
    public abstract String exon();

    @NotNull
    public abstract String grch37ChrStartRefAlt();

    @NotNull
    public abstract String nucleotideChange();

    @Nullable
    public abstract String aa();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract List<String> protCoords();
}
