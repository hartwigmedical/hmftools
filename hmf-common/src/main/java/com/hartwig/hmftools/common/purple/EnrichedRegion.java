package com.hartwig.hmftools.common.purple;

import java.io.Serializable;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class})
public abstract class EnrichedRegion implements GenomeRegion, Serializable {

    public abstract int mBAFCount();

    public abstract double mBAF();

    public abstract double tumorRatio();

    public abstract double normalRatio();
}
