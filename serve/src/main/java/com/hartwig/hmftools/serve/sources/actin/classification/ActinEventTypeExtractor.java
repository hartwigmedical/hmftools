package com.hartwig.hmftools.serve.sources.actin.classification;

import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;

import org.jetbrains.annotations.NotNull;

public class ActinEventTypeExtractor {

    @NotNull
    private static final EventClassifier CLASSIFIER = EventClassifierFactory.buildClassifier(ActinClassificationConfig.build());

    private ActinEventTypeExtractor() {
    }

    @NotNull
    public static EventType extractEventType(@NotNull ActinEntry entry) {
        String event = ActinEventExtractor.extractEvent(entry);

        return CLASSIFIER.determineType(entry.gene(), event);
    }
}
