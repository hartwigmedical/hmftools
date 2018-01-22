package com.hartwig.hmftools.common.cosmic.fusions;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class CosmicFusionModel {
    @NotNull
    public abstract List<CosmicFusionData> fusions();

    @NotNull
    public abstract List<CosmicPromiscuousFusionGene> promiscuousFivePrime();

    @NotNull
    public abstract List<CosmicPromiscuousFusionGene> promiscuousThreePrime();
}

