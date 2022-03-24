package com.hartwig.hmftools.common.cobalt;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ReadRatio extends GenomePosition
{

    double ratio();
}
