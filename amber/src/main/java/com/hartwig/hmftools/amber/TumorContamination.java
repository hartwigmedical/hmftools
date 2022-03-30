package com.hartwig.hmftools.amber;

import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface TumorContamination extends GenomePosition
{
    @NotNull
    BaseDepthData normal();

    @NotNull
    BaseDepthData tumor();
}
