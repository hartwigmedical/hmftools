package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

final class SignatureClassifier {

    private static final Set<String> SIGNATURES = Sets.newHashSet("Microsatellite Instability-High");

    private SignatureClassifier() {
    }

    public static boolean isSignature(@NotNull String event) {
        return SIGNATURES.contains(event);
    }
}
