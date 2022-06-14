package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class AmplificationMatcher implements EventMatcher {

    @NotNull
    private final Set<String> amplificationKeywords;
    @NotNull
    private final Set<String> amplificationKeyPhrases;


    AmplificationMatcher(@NotNull final Set<String> amplificationKeywords, @NotNull final Set<String> amplificationKeyPhrases) {
        this.amplificationKeywords = amplificationKeywords;
        this.amplificationKeyPhrases = amplificationKeyPhrases;
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

        return false;
    }
}
