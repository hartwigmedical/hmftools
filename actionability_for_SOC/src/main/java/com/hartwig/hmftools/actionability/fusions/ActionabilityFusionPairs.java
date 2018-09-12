package com.hartwig.hmftools.actionability.fusions;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class ActionabilityFusionPairs {

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
    abstract String levelSource();

    @NotNull
    abstract String hmfLevel();

    @NotNull
    abstract String evidenceType();

    @NotNull
    abstract String significanceSource();

    @NotNull
    abstract String hmfResponse();
}
