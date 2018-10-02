package com.hartwig.hmftools.actionability.cancerTypeMapping;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class CancerTypeMapping {

    @Nullable
    abstract String primaryTumorLocation();

    @Nullable
    abstract String doids();
}
