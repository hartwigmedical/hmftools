package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

public class GeneLevelMatcher implements EventMatcher {

    @NotNull
    private final Set<String> blackListKeyPhrases;
    @NotNull
    private final Set<String> genericGeneLevelKeyPhrases;
    @NotNull
    private final Set<String> activatingGeneLevelKeyPhrases;
    @NotNull
    private final Set<String> inactivatingGeneLevelKeyPhrases;

    public GeneLevelMatcher(@NotNull final Set<String> blackListKeyPhrases, @NotNull final Set<String> genericGeneLevelKeyPhrases,
            @NotNull final Set<String> activatingGeneLevelKeyPhrases, @NotNull final Set<String> inactivatingGeneLevelKeyPhrases) {
        this.blackListKeyPhrases = blackListKeyPhrases;
        this.genericGeneLevelKeyPhrases = genericGeneLevelKeyPhrases;
        this.activatingGeneLevelKeyPhrases = activatingGeneLevelKeyPhrases;
        this.inactivatingGeneLevelKeyPhrases = inactivatingGeneLevelKeyPhrases;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyPhrase : blackListKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return false;
            }
        }

        for (String keyPhrase : genericGeneLevelKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        for (String keyPhrase : activatingGeneLevelKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        for (String keyPhrase : inactivatingGeneLevelKeyPhrases) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        // If the event matches the gene we assume its a gene level event
        return event.trim().equals(gene);
    }
}
