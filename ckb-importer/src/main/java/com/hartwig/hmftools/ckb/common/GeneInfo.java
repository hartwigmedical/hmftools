package com.hartwig.hmftools.ckb.common;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneInfo {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract List<String> terms();
}
