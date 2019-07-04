package com.hartwig.hmftools.linx.visualiser.data;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionBuilder;

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

    public abstract String name();

    public abstract long namePosition();

}
