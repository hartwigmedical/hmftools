package com.hartwig.hmftools.common.actionability.fusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class ActionabilityFusionPairs {

    @NotNull
    abstract String fiveGene();

    @NotNull
    abstract String threeGene();

    @NotNull
    abstract String source();

    @NotNull
    abstract String reference();

    @NotNull
    abstract String drugsName();

    @NotNull
    abstract String drugsType();

    @NotNull
    abstract String cancerType();

    @NotNull
    abstract String hmfLevel();

    @NotNull
    abstract String hmfResponse();
}
