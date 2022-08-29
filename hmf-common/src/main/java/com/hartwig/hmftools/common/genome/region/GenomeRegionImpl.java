package com.hartwig.hmftools.common.genome.region;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true, passAnnotations = { NotNull.class, Nullable.class})
public interface GenomeRegionImpl extends GenomeRegion
{
}
