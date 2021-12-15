package com.hartwig.hmftools.serve.sources.actin.classification;

import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActinEventAndGeneExtractor {

    private ActinEventAndGeneExtractor() {
    }

    @NotNull
    public static String extractGene(@NotNull ActinEntry entry) {
        // TODO
        return Strings.EMPTY;
    }

    @NotNull
    public static String extractEvent(@NotNull ActinEntry entry) {
        // TODO
        switch (entry.rule()) {
            case ACTIVATION_OF_GENE_X:
                return Strings.EMPTY;
            case INACTIVATION_OF_GENE_X:
                return Strings.EMPTY;
            case ACTIVATING_MUTATION_IN_GENE_X:
                return Strings.EMPTY;
            case MUTATION_IN_GENE_X_OF_TYPE_Y:
                return Strings.EMPTY;
            case INACTIVATING_MUTATION_IN_GENE_X:
                return Strings.EMPTY;
            case AMPLIFICATION_OF_GENE_X:
                return Strings.EMPTY;
            case DELETION_OF_GENE_X:
                return Strings.EMPTY;
            case ACTIVATING_FUSION_IN_GENE_X:
                return Strings.EMPTY;
            case SPECIFIC_FUSION_X:
                return Strings.EMPTY;
            default: {
                throw new IllegalStateException("Unrecognized event: " + entry.rule());
            }
        }
    }
}
