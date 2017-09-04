package com.hartwig.hmftools.common.cosmic;

import java.util.Map;
import java.util.stream.Collectors;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class CosmicModel {
    @NotNull
    private final Map<String, CosmicData> geneData;
    private final Map<String, CosmicData> entrezMap;

    CosmicModel(@NotNull final Map<String, CosmicData> geneData) {
        this.geneData = geneData;
        entrezMap = geneData.entrySet().stream().collect(Collectors.toMap(v -> v.getValue().entrezId(), Map.Entry::getValue, (a, b) -> a));
    }

    @NotNull
    public Map<String, CosmicData> data() {
        return geneData;
    }

    @NotNull
    public Map<String, CosmicData> getEntrezMap() {
        return entrezMap;
    }

    @NotNull
    public String getRoleForGene(@NotNull final String gene) {
        final CosmicData geneCosmicData = geneData.get(gene);
        return geneCosmicData != null ? geneCosmicData.role() : Strings.EMPTY;
    }
}
