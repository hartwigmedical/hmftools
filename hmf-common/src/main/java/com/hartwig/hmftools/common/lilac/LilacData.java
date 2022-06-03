package com.hartwig.hmftools.common.lilac;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LilacData {

    @NotNull
    public abstract String qc();

    @NotNull
    public abstract List<LilacAllele> alleles();
}
