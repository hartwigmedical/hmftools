package com.hartwig.hmftools.linx.visualiser.data;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionBuilder;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Gene implements GenomeRegion
{
    public abstract static class Builder implements GenomeRegionBuilder<Gene>
    {
    }

    @NotNull
    public abstract String name();

    @NotNull
    public abstract VisGeneAnnotationType type();

    @NotNull
    public abstract String transcript();

    public abstract long namePosition();

    public abstract int strand();

}
