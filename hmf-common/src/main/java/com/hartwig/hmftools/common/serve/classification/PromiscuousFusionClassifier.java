package com.hartwig.hmftools.common.serve.classification;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public class PromiscuousFusionClassifier implements EventMatcher {

    private static final Set<String> FUSION_KEYWORDS =
            Sets.newHashSet("Fusion", "fusion", "FUSION", "Fusions", "FUSIONS", "REARRANGEMENT", "rearrange");

    @NotNull
    public static EventMatcher create(@NotNull List<EventMatcher> noMatchEventMatchers) {
        return new CompositeEventMatcher(noMatchEventMatchers, new PromiscuousFusionClassifier());
    }

    private PromiscuousFusionClassifier() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyword : FUSION_KEYWORDS) {
            if (event.contains(keyword) && !FusionPairClassifier.isFusionPair(event)) {
                return true;
            }
        }

        return false;
    }
}
