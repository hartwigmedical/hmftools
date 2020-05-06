package com.hartwig.hmftools.serve.actionability.variant;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class VariantEvidenceAnalyzer {

    @NotNull
    private final List<ActionableVariant> variant;

    VariantEvidenceAnalyzer(@NotNull final List<ActionableVariant> variant) {
        this.variant = variant;
    }
}
