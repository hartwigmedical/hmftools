package com.hartwig.hmftools.common.cosmic;

import java.util.Map;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class CosmicModel {

    @NotNull
    private final Map<String, CosmicData> geneData;

    CosmicModel(@NotNull final Map<String, CosmicData> geneData) {
        this.geneData = geneData;
    }

    @NotNull
    public Map<String, CosmicData> data() {
        return geneData;
    }

    @NotNull
    public String getRoleForGene(@NotNull final String gene) {
        final CosmicData geneCosmicData = geneData.get(gene);
        return geneCosmicData != null ? geneCosmicData.role() : Strings.EMPTY;
    }
}
