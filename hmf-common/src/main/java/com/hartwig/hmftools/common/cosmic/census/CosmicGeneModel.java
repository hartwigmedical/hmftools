package com.hartwig.hmftools.common.cosmic.census;

import java.util.Map;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class CosmicGeneModel {

    @NotNull
    private final Map<String, CosmicGeneData> geneData;

    CosmicGeneModel(@NotNull final Map<String, CosmicGeneData> geneData) {
        this.geneData = geneData;
    }

    @NotNull
    public Map<String, CosmicGeneData> data() {
        return geneData;
    }

    @NotNull
    public String getRoleForGene(@NotNull final String gene) {
        final CosmicGeneData geneCosmicGeneData = geneData.get(gene);
        return geneCosmicGeneData != null ? geneCosmicGeneData.role() : Strings.EMPTY;
    }
}
