package com.hartwig.hmftools.vicc.annotation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.ExclusiveEventClassifier;

import org.jetbrains.annotations.NotNull;

class DeletionClassifier implements EventClassifier {

    private static final Set<String> DELETION_KEYWORDS =
            Sets.newHashSet("Deletion", "deletion", "DELETION", "del", "undexpression", "UNDEREXPRESSION", "loss", "LOSS");

    private static final Set<String> DELETION_KEY_PHRASES = Sets.newHashSet("dec exp", "Copy Number Loss");

    private static final Set<String> KEYWORDS_TO_SKIP_FOR_DELETION = Sets.newHashSet("exon", "EXON", "Exon", "Ex19", "inframe");

    @NotNull
    public static EventClassifier create(@NotNull List<EventClassifier> excludingEventClassifiers) {
        return new ExclusiveEventClassifier(excludingEventClassifiers, new DeletionClassifier());
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
