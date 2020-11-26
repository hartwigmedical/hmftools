package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

class SignatureMatcher implements EventMatcher {

    private static final Set<String> SIGNATURES = Sets.newHashSet("Microsatellite Instability-High");

    SignatureMatcher() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        return SIGNATURES.contains(event);
    }
}
