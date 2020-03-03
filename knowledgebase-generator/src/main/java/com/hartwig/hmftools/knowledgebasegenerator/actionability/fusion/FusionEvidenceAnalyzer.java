package com.hartwig.hmftools.knowledgebasegenerator.actionability.fusion;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class FusionEvidenceAnalyzer {

    @NotNull
    private final List<ActionableFusion> fusions;

    FusionEvidenceAnalyzer(@NotNull final List<ActionableFusion> fusions) {
        this.fusions = fusions;
    }
}
