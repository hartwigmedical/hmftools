package com.hartwig.hmftools.ckb.classification;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.jetbrains.annotations.NotNull;

public final class EventTypeExtractor {

    @NotNull
    private static final EventClassifier CLASSIFIER = EventClassifierFactory.buildClassifier(CkbClassificationConfig.build());

    private EventTypeExtractor() {
    }

    @NotNull
    public static EventType classify(@NotNull CkbEntry entry) {
        if (entry.variants().size() > 1) {
            return EventType.COMBINED;
        } else {
            String event;
            Variant variant = entry.variants().get(0);
            if (variant.variant().equals("fusion") && variant.impact() != null && variant.impact().equals("fusion")) {
                event = "fusion promiscuous";
            } else if (variant.impact() != null && variant.impact().equals("fusion")) {
                event = variant.variant().replaceAll("\\s+", "") + " fusion";
            } else if (variant.variant().contains("exon")) {
                event = variant.variant().replace("exon", "exon ");
            } else {
                event = variant.variant();
            }

            return CLASSIFIER.determineType(variant.gene().geneSymbol(), event);
        }
    }
}
