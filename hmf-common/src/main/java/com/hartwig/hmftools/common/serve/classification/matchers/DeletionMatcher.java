package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class DeletionMatcher implements EventMatcher {

    @NotNull
    private final Set<String> blackListKeyPhrases;
    @NotNull
    private final Set<String> deletionKeywords;
    @NotNull
    private final Set<String> deletionKeyPhrases;

    public DeletionMatcher(@NotNull final Set<String> blackListKeyPhrases, @NotNull final Set<String> deletionKeywords,
            @NotNull final Set<String> deletionKeyPhrases) {
        this.blackListKeyPhrases = blackListKeyPhrases;
        this.deletionKeywords = deletionKeywords;
        this.deletionKeyPhrases = deletionKeyPhrases;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String blacklistPhrase : blackListKeyPhrases) {
            if (event.contains(blacklistPhrase)) {
                return false;
            }
        }

        String[] wordsDel = event.split(" ");
        for (String keyword : deletionKeywords) {
            for (String word : wordsDel) {
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
