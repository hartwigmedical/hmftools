package com.hartwig.hmftools.serve.sources.actin.classification;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActinEventAndGeneExtractor {

    private ActinEventAndGeneExtractor() {
    }

    @NotNull
    public static String extractGene(@NotNull ActinEntry entry) {
        return entry.parameters().get(0);
    }

    @NotNull
    public static List<String> extractEvent(@NotNull ActinEntry entry) {
        switch (entry.rule()) {
            case ACTIVATION_OF_GENE_X:
                return Lists.newArrayList("activation", "amplification");
            case INACTIVATION_OF_GENE_X:
                return Lists.newArrayList("inactivation", "deletion");
            case ACTIVATING_MUTATION_IN_GENE_X:
                return Lists.newArrayList("activation");
            case MUTATION_IN_GENE_X_OF_TYPE_Y:
                if (entry.parameters().size() == 2) {
                    return Lists.newArrayList(entry.parameters().get(1));
                } else {
                    return Lists.newArrayList(Strings.EMPTY);
                }
            case INACTIVATING_MUTATION_IN_GENE_X:
                return Lists.newArrayList("inactivation");
            case AMPLIFICATION_OF_GENE_X:
                return Lists.newArrayList("amplification");
            case DELETION_OF_GENE_X:
                return Lists.newArrayList("deletion");
            case ACTIVATING_FUSION_IN_GENE_X:
                return Lists.newArrayList("fusion");
            case SPECIFIC_FUSION_X:
                return Lists.newArrayList(entry.parameters().get(0).replace("_", "-"));
            default: {
                throw new IllegalStateException("Unrecognized event: " + entry.rule());
            }
        }
    }
}
