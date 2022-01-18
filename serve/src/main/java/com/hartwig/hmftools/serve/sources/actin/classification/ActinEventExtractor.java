package com.hartwig.hmftools.serve.sources.actin.classification;

import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.jetbrains.annotations.NotNull;

public final class ActinEventExtractor {

    private ActinEventExtractor() {
    }

    @NotNull
    public static String extractEvent(@NotNull ActinEntry entry) {
        switch (entry.rule()) {
            case ACTIVATING_MUTATION_IN_GENE_X:
            case ACTIVATION_OF_GENE_X:
                return ActinKeywords.ACTIVATION;
            case INACTIVATING_MUTATION_IN_GENE_X:
            case INACTIVATION_OF_GENE_X:
                return ActinKeywords.INACTIVATION;
            case MUTATION_IN_GENE_X_OF_TYPE_Y: {
                String mutation = entry.mutation();
                if (mutation == null) {
                    throw new IllegalStateException("No mutation provided in ACTIN entry: " + entry);
                }
                return mutation;
            }
            case AMPLIFICATION_OF_GENE_X:
                return ActinKeywords.AMPLIFICATION;
            case DELETION_OF_GENE_X:
                return ActinKeywords.DELETION;
            case ACTIVATING_FUSION_IN_GENE_X:
                return ActinKeywords.PROMISCUOUS_FUSION;
            case SPECIFIC_FUSION_X:
                return entry.gene() + " fusion";
            default: {
                throw new IllegalStateException("Unrecognized event: " + entry.rule());
            }
        }
    }
}
