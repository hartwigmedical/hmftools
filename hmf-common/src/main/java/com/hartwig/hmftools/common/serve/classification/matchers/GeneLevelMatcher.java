package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

public class GeneLevelMatcher implements EventMatcher {

    @NotNull
    private final String exonKeyword;
    @NotNull
    private final Set<String> genericGeneLevelKeywords;
    @NotNull
    private final Set<String> activatingGeneLevelKeywords;
    @NotNull
    private final Set<String> inactivatingGeneLevelKeywords;

    public GeneLevelMatcher(@NotNull final String exonKeyword, @NotNull final Set<String> genericGeneLevelKeywords,
            @NotNull final Set<String> activatingGeneLevelKeywords, @NotNull final Set<String> inactivatingGeneLevelKeywords) {
        this.exonKeyword = exonKeyword;
        this.genericGeneLevelKeywords = genericGeneLevelKeywords;
        this.activatingGeneLevelKeywords = activatingGeneLevelKeywords;
        this.inactivatingGeneLevelKeywords = inactivatingGeneLevelKeywords;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        if (event.toLowerCase().contains(exonKeyword)) {
            return false;
        }

        for (String keyword : genericGeneLevelKeywords) {
            if (event.contains(keyword)) {
                return true;
            }
        }

        for (String keyword : activatingGeneLevelKeywords) {
            if (event.contains(keyword)) {
                return true;
            }
        }

        for (String keyword : inactivatingGeneLevelKeywords) {
            if (event.contains(keyword)) {
                return true;
            }
        }

        return event.trim().equals(gene);
    }
}
