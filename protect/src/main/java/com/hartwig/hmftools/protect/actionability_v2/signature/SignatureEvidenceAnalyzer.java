package com.hartwig.hmftools.protect.actionability_v2.signature;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class SignatureEvidenceAnalyzer {

    @NotNull
    private final List<ActionableSignature> signature;

    SignatureEvidenceAnalyzer(@NotNull final List<ActionableSignature> signature) {
        this.signature = signature;
    }
}
