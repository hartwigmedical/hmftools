package com.hartwig.hmftools.linx.visualiser.data;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionBuilder;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ProteinDomain implements GenomeRegion
{
    public abstract static class Builder implements GenomeRegionBuilder<ProteinDomain>
    {
    }

    @NotNull
    public abstract String sampleId();

    public abstract int clusterId();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String transcript();
}
