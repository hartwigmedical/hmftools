package com.hartwig.hmftools.common.serve.classification;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

class DeletionClassifier implements EventMatcher {

    private static final Set<String> DELETION_KEYWORDS =
            Sets.newHashSet("Deletion", "deletion", "DELETION", "del", "undexpression", "UNDEREXPRESSION", "loss", "LOSS");

    private static final Set<String> DELETION_KEY_PHRASES = Sets.newHashSet("dec exp", "Copy Number Loss");

    private static final Set<String> KEYWORDS_TO_SKIP_FOR_DELETION = Sets.newHashSet("exon", "EXON", "Exon", "Ex19", "inframe");

    @NotNull
    public static EventMatcher create(@NotNull List<EventMatcher> noMatchEventMatchers) {
        return new CompositeEventMatcher(noMatchEventMatchers, new DeletionClassifier());
    }

    private DeletionClassifier() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String skipTerm : KEYWORDS_TO_SKIP_FOR_DELETION) {
            if (event.contains(skipTerm)) {
                return false;
            }
        }

        String[] words = event.split(" ");
        for (String keyword : DELETION_KEYWORDS) {
            for (String word : words) {
                if (word.equals(keyword)) {
                    return true;
                }
            }
        }

        for (String keyPhrase : DELETION_KEY_PHRASES) {
            if (event.contains(keyPhrase)) {
                return true;
            }
        }

        return false;
    }
}
