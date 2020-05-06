package com.hartwig.hmftools.knowledgebasegenerator.actionability.range;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class RangeEvidenceAnalyzer {

    @NotNull
    private final List<ActionableRange> range;

    RangeEvidenceAnalyzer(@NotNull final List<ActionableRange> range) {
        this.range = range;
    }
}
