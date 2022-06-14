package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class AmplificationMatcher implements EventMatcher {

    @NotNull
    private final Set<String> amplificationKeywords;
    @NotNull
    private final Set<String> amplificationKeyPhrases;
    @NotNull
    private final Set<String> overexpressionKeywords;
    @NotNull
    private final Set<String> overexpressionKeyPhrases;

    AmplificationMatcher(@NotNull final Set<String> amplificationKeywords, @NotNull final Set<String> amplificationKeyPhrases,
            @NotNull final Set<String> overexpressionKeywords, @NotNull final Set<String> overexpressionKeyPhrases) {
        this.amplificationKeywords = amplificationKeywords;
        this.amplificationKeyPhrases = amplificationKeyPhrases;
        this.overexpressionKeywords = overexpressionKeywords;
        this.overexpressionKeyPhrases = overexpressionKeyPhrases;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        String[] wordsAmp = event.split(" ");
        for (String keyword : amplificationKeywords) {
            for (String word : wordsAmp) {
                if (word.equals(keyword)) {
                    return true;
                }
            }
        }

        for (String keyPhrase : amplificationKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

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
