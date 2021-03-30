package com.hartwig.hmftools.ckb.classification;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
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
    public static EventType classify(@NotNull CkbEntry entry) {

        String gene = entry.variants().get(0).gene().geneSymbol();

        String profileName;
        EventType type;
        if (entry.variants().size() > 1) {
            return EventType.COMBINED;
        } else {
            if (entry.variants().get(0).variant().equals("fusion") && entry.variants().get(0).impact() != null && entry.variants()
                    .get(0)
                    .impact()
                    .equals("fusion")) {
                profileName = "fusion promisuous";
            } else if (entry.variants().get(0).impact() != null && entry.variants().get(0).impact().equals("fusion")) {
                profileName = entry.variants().get(0).variant().replaceAll("\\s+","") + " fusion";
            } else if (entry.variants().get(0).variant().contains("exon")) {
                profileName = entry.variants().get(0).variant().replace("exon", "exon ");
            }
            else {
                profileName = entry.variants().get(0).variant() ;
            }

            return CLASSIFIER.determineType(gene, profileName);
        }

    }
}
