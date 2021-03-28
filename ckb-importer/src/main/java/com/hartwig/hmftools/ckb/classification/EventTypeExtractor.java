package com.hartwig.hmftools.ckb.classification;

import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.jetbrains.annotations.NotNull;

public class EventTypeExtractor {

    @NotNull
    private static final EventClassifier CLASSIFIER = EventClassifierFactory.buildClassifier(CkbClassificationConfig.build());

    private EventTypeExtractor() {
    }

    @NotNull
    public static EventType classify(@NotNull String gene, @NotNull String variant) {
        return CLASSIFIER.determineType(gene, variant);
    }
}
