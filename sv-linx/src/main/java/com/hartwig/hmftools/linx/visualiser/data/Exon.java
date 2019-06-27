package com.hartwig.hmftools.linx.visualiser.data;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionBuilderI;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Exon implements GenomeRegion
{
    public abstract static class Builder implements GenomeRegionBuilderI<Exon>
    {
    }

    public abstract String sampleId();

    public abstract int clusterId();

    public abstract String gene();

    public abstract int rank();

}
