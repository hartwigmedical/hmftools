package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class SignatureMatcher implements EventMatcher {

    @NotNull
    private final Set<String> signatureEvents;

    SignatureMatcher(@NotNull final Set<String> signatureEvents) {
        this.signatureEvents = signatureEvents;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        return signatureEvents.contains(event);
    }
}
