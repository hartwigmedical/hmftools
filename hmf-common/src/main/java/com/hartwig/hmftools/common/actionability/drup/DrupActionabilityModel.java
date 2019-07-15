package com.hartwig.hmftools.common.actionability.drup;

import java.util.Set;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrupActionabilityModel {

    @NotNull
    public abstract Set<String> actionableGenes();

}
