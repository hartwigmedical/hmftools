package com.hartwig.hmftools.common.serve.classification;

import com.hartwig.hmftools.common.serve.classification.matchers.MatcherFactory;

import org.jetbrains.annotations.NotNull;

public final class EventClassifierFactory {

    private EventClassifierFactory() {
    }

    @NotNull
    public static EventClassifier buildClassifier(@NotNull EventPreprocessor proteinAnnotationPreprocessor) {
        return new EventClassifier(MatcherFactory.buildMatcherMap(proteinAnnotationPreprocessor));
    }
}
