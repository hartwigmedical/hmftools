package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class DeletionMatcher implements EventMatcher {

    @NotNull
    private final Set<String> deletionKeywords;
    @NotNull
    private final Set<String> deletionKeyPhrases;
    @NotNull
    private final Set<String> deletionKeyPhrasesToSkip;

    DeletionMatcher(@NotNull final Set<String> deletionKeywords, @NotNull final Set<String> deletionKeyPhrases,
            @NotNull final Set<String> deletionKeyPhrasesToSkip) {
        this.deletionKeywords = deletionKeywords;
        this.deletionKeyPhrases = deletionKeyPhrases;
        this.deletionKeyPhrasesToSkip = deletionKeyPhrasesToSkip;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String skipTerm : deletionKeyPhrasesToSkip) {
            if (event.contains(skipTerm)) {
                return false;
            }
        }

        String[] words = event.split(" ");
        for (String keyword : deletionKeywords) {
            for (String word : words) {
                if (word.equals(keyword)) {
                    return true;
                }
            }
        }

        for (String keyPhrase : deletionKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        return false;
    }
}
