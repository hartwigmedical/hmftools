package com.hartwig.hmftools.common.amber;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface TumorContamination extends GenomePosition {

    @NotNull
    BaseDepth normal();

    @NotNull
    BaseDepth tumor();
}
