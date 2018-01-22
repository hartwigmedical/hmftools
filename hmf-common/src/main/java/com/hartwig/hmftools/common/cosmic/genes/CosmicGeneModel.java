package com.hartwig.hmftools.common.cosmic.genes;

import java.util.Map;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class CosmicGeneModel {

    @NotNull
    public abstract Map<String, CosmicGeneData> data();

    @NotNull
    public String getRoleForGene(@NotNull final String gene) {
        final CosmicGeneData cosmicGeneData = data().get(gene);
        return cosmicGeneData != null ? cosmicGeneData.role() : Strings.EMPTY;
    }
}
