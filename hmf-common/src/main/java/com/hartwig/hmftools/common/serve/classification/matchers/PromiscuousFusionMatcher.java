package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;

import org.jetbrains.annotations.NotNull;

class PromiscuousFusionMatcher implements EventMatcher {

    private static final Set<String> FUSION_KEYWORDS =
            Sets.newHashSet("Fusion", "fusion", "FUSION", "Fusions", "FUSIONS", "REARRANGEMENT", "rearrange");

    @NotNull
    public static EventMatcher create(@NotNull List<EventMatcher> noMatchEventMatchers) {
        return new CompositeEventMatcher(noMatchEventMatchers, new PromiscuousFusionMatcher());
    }

    @VisibleForTesting
    PromiscuousFusionMatcher() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyword : FUSION_KEYWORDS) {
            if (event.contains(keyword) && !FusionPairMatcher.isFusionPair(event)) {
                return true;
            }
        }

        return false;
    }
}
