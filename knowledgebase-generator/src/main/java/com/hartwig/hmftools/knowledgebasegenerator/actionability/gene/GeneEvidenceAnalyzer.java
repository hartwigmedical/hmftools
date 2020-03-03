package com.hartwig.hmftools.knowledgebasegenerator.actionability.gene;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class GeneEvidenceAnalyzer {

    @NotNull
    private final List<ActionableGene> gene;

    GeneEvidenceAnalyzer(@NotNull final List<ActionableGene> gene) {
        this.gene = gene;
    }
}
