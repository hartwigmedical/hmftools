package com.hartwig.hmftools.linx.visualiser.data;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionBuilder;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Exon implements GenomeRegion
{
    public abstract static class Builder implements GenomeRegionBuilder<Exon>
    {
    }

    @NotNull
    public abstract String sampleId();

    public abstract int clusterId();

    public abstract VisGeneAnnotationType type();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcript();

    public abstract int rank();

}
