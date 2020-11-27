package com.hartwig.hmftools.common.serve.classification;

import com.hartwig.hmftools.common.serve.classification.matchers.EventMatcherFactory;

import org.jetbrains.annotations.NotNull;

public final class EventClassifierFactory {

    private EventClassifierFactory() {
    }

    @NotNull
    public static EventClassifier buildClassifier(@NotNull EventClassifierConfig config) {
        return new EventClassifier(EventMatcherFactory.buildMatcherMap(config));
    }
}
