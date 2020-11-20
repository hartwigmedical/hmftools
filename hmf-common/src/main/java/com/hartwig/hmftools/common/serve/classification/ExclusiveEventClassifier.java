package com.hartwig.hmftools.common.serve.classification;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class ExclusiveEventClassifier implements EventClassifier {

    @NotNull
    private final List<EventClassifier> excludingEventClassifiers;
    @NotNull
    private final EventClassifier inclusiveEventClassifier;

    public ExclusiveEventClassifier(@NotNull final List<EventClassifier> excludingEventClassifiers,
            @NotNull final EventClassifier inclusiveEventClassifier) {
        this.excludingEventClassifiers = excludingEventClassifiers;
        this.inclusiveEventClassifier = inclusiveEventClassifier;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        for (EventClassifier classifier : excludingEventClassifiers) {
            if (classifier.matches(gene, event)) {
                return false;
            }
        }
        return inclusiveEventClassifier.matches(gene, event);
    }
}
