package com.hartwig.hmftools.protect.actionability_v2.range;

import java.util.List;


import org.jetbrains.annotations.NotNull;

public class RangeEvidenceAnalyzer {

    @NotNull
    private final List<ActionableRange> range;

    RangeEvidenceAnalyzer(@NotNull final List<ActionableRange> range) {
        this.range = range;
    }
}
