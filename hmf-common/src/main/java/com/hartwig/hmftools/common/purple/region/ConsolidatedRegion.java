package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ConsolidatedRegion implements GenomeRegion {

    public abstract double averageBAF();

    public abstract double averageRatioOfRatios();
}
