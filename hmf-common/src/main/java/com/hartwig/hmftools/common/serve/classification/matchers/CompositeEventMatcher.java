package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class CompositeEventMatcher implements EventMatcher {

    @NotNull
    private final List<EventMatcher> noMatchEventMatchers;
    @NotNull
    private final EventMatcher finalEventMatcher;

    CompositeEventMatcher(@NotNull final List<EventMatcher> noMatchEventMatchers, @NotNull final EventMatcher finalEventMatcher) {
        this.noMatchEventMatchers = noMatchEventMatchers;
        this.finalEventMatcher = finalEventMatcher;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (EventMatcher classifier : noMatchEventMatchers) {
            if (classifier.matches(gene, event)) {
                return false;
            }
        }
        return finalEventMatcher.matches(gene, event);
    }
}
