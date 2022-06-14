package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class OverExpressionMatcher implements EventMatcher {

    @NotNull
    private final Set<String> overexpressionKeywords;
    @NotNull
    private final Set<String> overexpressionKeyPhrases;

    OverExpressionMatcher(@NotNull final Set<String> overexpressionKeywords, @NotNull final Set<String> overexpressionKeyPhrases) {
        this.overexpressionKeywords = overexpressionKeywords;
        this.overexpressionKeyPhrases = overexpressionKeyPhrases;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {

        String[] wordsOver = event.split(" ");
        for (String keyword : overexpressionKeywords) {
            for (String word : wordsOver) {
                if (word.equals(keyword)) {
                    return true;
                }
            }
        }

        for (String keyPhrase : overexpressionKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        return false;
    }
}
