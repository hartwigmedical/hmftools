package com.hartwig.hmftools.common.cosmic.fusions;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class CosmicFusionData {
    @NotNull
    public abstract String fiveGene();

    @Nullable
    public abstract String fiveTranscript();

    @NotNull
    public abstract String threeGene();

    @Nullable
    public abstract String threeTranscript();

    @NotNull
    public abstract String cosmicURL();
}
