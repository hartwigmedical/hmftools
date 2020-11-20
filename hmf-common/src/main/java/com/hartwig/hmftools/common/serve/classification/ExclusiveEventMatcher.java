package com.hartwig.hmftools.common.serve.classification;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class ExclusiveEventMatcher implements EventMatcher {

    @NotNull
    private final List<EventMatcher> excludingEventMatchers;
    @NotNull
    private final EventMatcher inclusiveEventMatcher;

    public ExclusiveEventMatcher(@NotNull final List<EventMatcher> excludingEventMatchers,
            @NotNull final EventMatcher inclusiveEventMatcher) {
        this.excludingEventMatchers = excludingEventMatchers;
        this.inclusiveEventMatcher = inclusiveEventMatcher;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (EventMatcher classifier : excludingEventMatchers) {
            if (classifier.matches(gene, event)) {
                return false;
            }
        }
        return inclusiveEventMatcher.matches(gene, event);
    }
}
