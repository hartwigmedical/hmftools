package com.hartwig.hmftools.common.serve.classification;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

class GeneRangeExonMatcher implements EventMatcher {

    private static final String EXON_KEYWORD = "exon";
    private static final Set<String> EXON_RANGE_EXACT_TERMS = Sets.newHashSet("RARE EX 18-21 MUT");

    private static final Set<String> EXON_RANGE_KEYWORDS =
            Sets.newHashSet("deletion", "insertion", "proximal", "mutation", "splice site insertion", "frameshift");

    @NotNull
    public static EventMatcher create(@NotNull List<EventMatcher> noMatchEventMatchers) {
        return new CompositeEventMatcher(noMatchEventMatchers, new GeneRangeExonMatcher());
    }

    @VisibleForTesting
    GeneRangeExonMatcher() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        if (EXON_RANGE_EXACT_TERMS.contains(event)) {
            return true;
        } else {
            String lowerCaseEvent = event.toLowerCase();
            if (lowerCaseEvent.contains(EXON_KEYWORD)) {
                for (String keyword : EXON_RANGE_KEYWORDS) {
                    if (lowerCaseEvent.contains(keyword)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }
}
