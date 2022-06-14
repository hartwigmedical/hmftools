package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class UnderExpressionMatcher implements EventMatcher {

    @NotNull
    private final Set<String> underexpressionKeywords;
    @NotNull
    private final Set<String> underexpressionKeyPhrases;

    public UnderExpressionMatcher(@NotNull final Set<String> underexpressionKeywords,
            @NotNull final Set<String> underexpressionKeyPhrases) {
        this.underexpressionKeywords = underexpressionKeywords;
        this.underexpressionKeyPhrases = underexpressionKeyPhrases;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {

        String[] wordsUnder = event.split(" ");
        for (String keyword : underexpressionKeywords) {
            for (String word : wordsUnder) {
                if (word.equals(keyword)) {
                    return true;
                }
            }
        }

        for (String keyPhrase : underexpressionKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        return false;
    }
}
