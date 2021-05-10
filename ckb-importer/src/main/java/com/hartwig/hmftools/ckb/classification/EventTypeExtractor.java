package com.hartwig.hmftools.ckb.classification;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class EventTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(EventTypeExtractor.class);
    private static final EventClassifier CLASSIFIER = EventClassifierFactory.buildClassifier(CkbClassificationConfig.build());
    private static final EventAndGeneExtractor EXTRACTOR = new EventAndGeneExtractor();

    private EventTypeExtractor() {
    }

    @NotNull
    public static EventType classify(@NotNull CkbEntry entry) {
        int variantCount = entry.variants().size();
        if (variantCount > 1) {
            return EventType.COMBINED;
        } else if (variantCount == 1) {
            Variant variant = entry.variants().get(0);

            String gene = EXTRACTOR.extractGene(variant);
            String event = EXTRACTOR.extractEvent(variant);

            return CLASSIFIER.determineType(gene, event);
        } else {
            LOGGER.warn("CKB entry found with no variants: {}", entry);
            return EventType.UNKNOWN;
        }
    }
}
