package com.hartwig.hmftools.protect.actionability_v2.variant;

import java.util.List;


import org.jetbrains.annotations.NotNull;

public class VariantEvidenceAnalyzer {

    @NotNull
    private final List<ActionableVariant> variant;

    VariantEvidenceAnalyzer(@NotNull final List<ActionableVariant> variant) {
        this.variant = variant;
    }
}
