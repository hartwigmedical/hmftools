package com.hartwig.hmftools.knowledgebasegenerator.actionability.signature;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class SignatureEvidenceAnalyzer {

    @NotNull
    private final List<ActionableSignature> signature;

    SignatureEvidenceAnalyzer(@NotNull final List<ActionableSignature> signature) {
        this.signature = signature;
    }
}
