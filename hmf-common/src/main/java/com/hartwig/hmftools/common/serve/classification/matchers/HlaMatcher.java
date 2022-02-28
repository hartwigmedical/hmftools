package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.Set;

import org.jetbrains.annotations.NotNull;

public class HlaMatcher implements EventMatcher{

    @NotNull
    private final Set<String> hlaEvents;

    HlaMatcher(@NotNull final Set<String> hlaEvents) {
        this.hlaEvents = hlaEvents;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        return hlaEvents.contains(event);
    }
}
