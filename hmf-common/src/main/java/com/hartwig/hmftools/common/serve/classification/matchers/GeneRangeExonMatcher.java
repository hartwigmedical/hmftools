package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

class GeneRangeExonMatcher implements EventMatcher {

    @NotNull
    private final String exonKeyword;
    @NotNull
    private final Set<String> exonRangeEvents;
    @NotNull
    private final Set<String> exonRangeKeywords;

    GeneRangeExonMatcher(@NotNull final String exonKeyword, @NotNull final Set<String> exonRangeEvents,
            @NotNull final Set<String> exonRangeKeywords) {
        this.exonKeyword = exonKeyword;
        this.exonRangeEvents = exonRangeEvents;
        this.exonRangeKeywords = exonRangeKeywords;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        if (exonRangeEvents.contains(event)) {
            return true;
        } else {
            String lowerCaseEvent = event.toLowerCase();
            if (lowerCaseEvent.contains(exonKeyword)) {
                for (String keyword : exonRangeKeywords) {
                    if (lowerCaseEvent.contains(keyword)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }
}
