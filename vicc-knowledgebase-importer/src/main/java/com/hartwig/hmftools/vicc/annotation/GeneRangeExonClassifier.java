package com.hartwig.hmftools.vicc.annotation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.ExclusiveEventClassifier;

import org.jetbrains.annotations.NotNull;

class GeneRangeExonClassifier implements EventClassifier {

    private static final String EXON_KEYWORD = "exon";
    private static final Set<String> EXON_RANGE_EXACT_TERMS = Sets.newHashSet("RARE EX 18-21 MUT");

    private static final Set<String> EXON_RANGE_KEYWORDS =
            Sets.newHashSet("deletion", "insertion", "proximal", "mutation", "splice site insertion", "frameshift");

    @NotNull
    public static EventClassifier create(@NotNull List<EventClassifier> excludingEventClassifiers) {
        return new ExclusiveEventClassifier(excludingEventClassifiers, new GeneRangeExonClassifier());
    }

    private GeneRangeExonClassifier() {
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
