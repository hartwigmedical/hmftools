package com.hartwig.hmftools.vicc.datamodel.oncokb;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OncoKbGene {

    @NotNull
    public abstract String hugoSymbol();

    @NotNull
    public abstract List<String> geneAliases();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String entrezGeneId();

    @Nullable
    public abstract String curatedIsoform();

    @Nullable
    public abstract String curatedRefSeq();

    @NotNull
    public abstract String oncogene();

    @NotNull
    public abstract String tsg();

}
