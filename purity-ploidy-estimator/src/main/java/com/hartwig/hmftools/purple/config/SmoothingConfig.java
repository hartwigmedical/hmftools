package com.hartwig.hmftools.purple.config;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SmoothingConfig {

    int minDiploidTumorRatioCount();

    int minDiploidTumorRatioCountAtCentromere();
}
