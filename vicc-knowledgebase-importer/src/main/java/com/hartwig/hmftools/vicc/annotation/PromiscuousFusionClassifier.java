package com.hartwig.hmftools.vicc.annotation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.ExclusiveEventClassifier;

import org.jetbrains.annotations.NotNull;

public class PromiscuousFusionClassifier implements EventClassifier {

    private static final Set<String> FUSION_KEYWORDS =
            Sets.newHashSet("Fusion", "fusion", "FUSION", "Fusions", "FUSIONS", "REARRANGEMENT", "rearrange");

    @NotNull
    public static EventClassifier create(@NotNull List<EventClassifier> excludingEventClassifiers) {
        return new ExclusiveEventClassifier(excludingEventClassifiers, new PromiscuousFusionClassifier());
    }

    private PromiscuousFusionClassifier() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (String keyword : FUSION_KEYWORDS) {
            if (event.contains(keyword)) {
                return true;
            }
        }

        return false;
    }
}
