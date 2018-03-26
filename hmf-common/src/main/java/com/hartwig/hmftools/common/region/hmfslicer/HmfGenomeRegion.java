package com.hartwig.hmftools.common.region.hmfslicer;

import java.util.List;

import com.hartwig.hmftools.common.gene.TranscriptRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class HmfGenomeRegion implements TranscriptRegion {

    @NotNull
    public abstract String geneID();

    public abstract long geneStart();

    public abstract long geneEnd();

    public abstract long codingStart();

    public abstract long codingEnd();

    @NotNull
    public abstract Strand strand();

    @NotNull
    public abstract List<Integer> entrezId();

    @NotNull
    public abstract List<HmfExonRegion> exome();
}
