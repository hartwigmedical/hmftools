package com.hartwig.hmftools.linx.visualiser.data;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionBuilder;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class FusedExon implements GenomeRegion
{
    public abstract static class Builder implements GenomeRegionBuilder<FusedExon> {}

    public abstract String sampleId();
    public abstract int clusterId();
    public abstract String fusion();
    public abstract String gene();
    public abstract boolean isUpstream();
    public abstract String transcript();
    public abstract int unadjustedGeneStart();
    public abstract int geneStart();
    public abstract int geneEnd();
    public abstract int rank();
    public abstract boolean skipped();

}
