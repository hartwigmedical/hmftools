package com.hartwig.hmftools.serve.sources.actin.classification;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.jetbrains.annotations.NotNull;

public final class ActinEventExtractor {

    private ActinEventExtractor() {
    }

    @NotNull
    public static Set<String> extractEvents(@NotNull ActinEntry entry) {
        switch (entry.rule()) {
            case ACTIVATION_OR_AMPLIFICATION_OF_GENE_X: {
                return Sets.newHashSet(ActinKeywords.ACTIVATION, ActinKeywords.AMPLIFICATION);
            }
            case ACTIVATING_MUTATION_IN_GENE_X: {
                return Sets.newHashSet(ActinKeywords.ACTIVATION);
            }
            case INACTIVATION_OF_GENE_X: {
                return Sets.newHashSet(ActinKeywords.INACTIVATION);
            }
            case MSI_SIGNATURE: {
                return Sets.newHashSet(ActinKeywords.MSI_SIGNATURE);
            }
            case HRD_SIGNATURE: {
                return Sets.newHashSet(ActinKeywords.HRD_SIGNATURE);
            }
            case TMB_OF_AT_LEAST_X:
            case TML_OF_AT_LEAST_X:
            case TML_OF_AT_MOST_X:
            case HAS_HLA_A_TYPE_X:
            case MUTATION_IN_GENE_X_OF_TYPE_Y: {
                String mutation = entry.mutation();
                if (mutation == null) {
                    throw new IllegalStateException("No mutation provided in ACTIN entry: " + entry);
                }
                return Sets.newHashSet(mutation);
            }
            case AMPLIFICATION_OF_GENE_X: {
                return Sets.newHashSet(ActinKeywords.AMPLIFICATION);
            }
            case DELETION_OF_GENE_X: {
                return Sets.newHashSet(ActinKeywords.DELETION);
            }
            case FUSION_IN_GENE_X: {
                return Sets.newHashSet(ActinKeywords.PROMISCUOUS_FUSION);
            }
            case SPECIFIC_FUSION_OF_X_TO_Y: {
                String mutation = entry.mutation();
                if (mutation == null) {
                    throw new IllegalStateException("No mutation provided in ACTIN entry: " + entry);
                }
                return Sets.newHashSet(entry.mutation() + " fusion");
            }
            case WILDTYPE_OF_GENE_X: {
                return Sets.newHashSet(ActinKeywords.WILDTYPE);
            }
            default: {
                throw new IllegalStateException("Unrecognized event: " + entry.rule());
            }
        }
    }
}
