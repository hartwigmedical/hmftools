package com.hartwig.hmftools.common.region.hmfslicer;

import java.util.List;

import com.hartwig.hmftools.common.gene.GeneRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class HmfGenomeRegion implements GeneRegion {

    @NotNull
    public abstract String geneID();

    public abstract long geneStart();

    public abstract long geneEnd();

    @NotNull
    public abstract String chromosomeBand();

    @NotNull
    public abstract List<Integer> entrezId();

    @NotNull
    public abstract List<HmfExonRegion> exome();
}
