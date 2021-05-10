package com.hartwig.hmftools.vicc.annotation;

import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ViccEventTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(ViccEventTypeExtractor.class);

    private static final EventClassifier CLASSIFIER = EventClassifierFactory.buildClassifier(ViccClassificationConfig.build());

    private ViccEventTypeExtractor() {
    }

    @NotNull
    public static EventType extractType(@NotNull Feature feature) {
        String gene = feature.geneSymbol();
        if (gene == null) {
            LOGGER.debug("Skipping extraction for '{}' since gene is missing", feature.name());
            return EventType.UNKNOWN;
        } else {
            return CLASSIFIER.determineType(gene, feature.name());
        }
    }
}
