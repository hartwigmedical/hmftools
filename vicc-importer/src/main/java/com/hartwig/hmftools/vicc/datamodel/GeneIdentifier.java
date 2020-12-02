package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneIdentifier {

    @NotNull
    public abstract String symbol();

    @NotNull
    public abstract String entrezId();

    @Nullable
    public abstract String ensemblGeneId();
}
