package com.hartwig.hmftools.linx.visualiser.data;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionBuilderI;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CopyNumberAlteration implements GenomeRegion
{
    public abstract static class Builder implements GenomeRegionBuilderI<CopyNumberAlteration>
    {
    }

    public abstract String sampleId();

    public abstract double copyNumber();

    public abstract double baf();

    public double minorAllelePloidy()
    {
        return Math.max(0, (1 - baf()) * copyNumber());
    }
}
