package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

class AmplificationMatcher implements EventMatcher {

    private static final Set<String> AMPLIFICATION_KEYWORDS =
            Sets.newHashSet("Amplification", "amplification", "AMPLIFICATION", "amp", "overexpression", "OVEREXPRESSION", "Overexpression");

    private static final Set<String> AMPLIFICATION_KEY_PHRASES = Sets.newHashSet("over exp");

    AmplificationMatcher() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        String[] words = event.split(" ");
        for (String keyword : AMPLIFICATION_KEYWORDS) {
            for (String word : words) {
                if (word.equals(keyword)) {
                    return true;
                }
            }
        }

        for (String keyPhrase : AMPLIFICATION_KEY_PHRASES) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        return false;
    }
}
