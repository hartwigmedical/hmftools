package com.hartwig.hmftools.serve.sources.actin.classification;

import com.hartwig.hmftools.common.genome.genepanel.GeneNameMapping37to38;
import com.hartwig.hmftools.serve.sources.actin.ActinTrial;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ActinEventTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ActinEventTypeExtractor.class);

    @NotNull
    private final GeneNameMapping37to38 geneNameMapping;

    public ActinEventTypeExtractor() {
        this.geneNameMapping = GeneNameMapping37to38.loadFromEmbeddedResource();
    }

    @NotNull
    public String extractGene(@NotNull ActinTrial trial) {
        return Strings.EMPTY;
    }

    @NotNull
    public String extractEvent(@NotNull ActinTrial trial) {
        switch (trial.rule()) {
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
            case OVEREXPRESSION_OF_GENE_X:
                return Strings.EMPTY;
            case WILDTYPE_OF_GENE_X:
                return Strings.EMPTY;
            case MSI_SIGNATURE:
                return Strings.EMPTY;
            case HRD_SIGNATURE:
                return Strings.EMPTY;
            case TMB_OF_AT_LEAST_X:
                return Strings.EMPTY;
            case TML_OF_AT_LEAST_X:
                return Strings.EMPTY;
            case TML_OF_AT_MOST_X:
                return Strings.EMPTY;
            default: {
                throw new IllegalStateException("Unrecognized event: " + trial.rule());
            }
        }
    }
}
