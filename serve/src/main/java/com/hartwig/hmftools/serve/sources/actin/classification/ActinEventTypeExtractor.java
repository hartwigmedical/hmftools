package com.hartwig.hmftools.serve.sources.actin.classification;

import java.util.List;

import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;
import com.hartwig.hmftools.serve.sources.ckb.CkbExtractor;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ActinEventTypeExtractor {

    @NotNull
    private static final EventClassifier CLASSIFIER = EventClassifierFactory.buildClassifier(ActinClassificationConfig.build());
    private static final Logger LOGGER = LogManager.getLogger(ActinEventTypeExtractor.class);

    private ActinEventTypeExtractor() {
    }

    @NotNull
    public static List<EventType> classify(@NotNull ActinEntry actinEntry) {
        List<EventType> eventType = Lists.newArrayList();
        String gene = ActinEventAndGeneExtractor.extractGene(actinEntry);
        List<String> events = ActinEventAndGeneExtractor.extractEvent(actinEntry);

        if (actinEntry.rule() == ActinRule.ACTIVATION_OF_GENE_X || actinEntry.rule() == ActinRule.INACTIVATION_OF_GENE_X) {
            for (String event : events) {
                eventType.add(CLASSIFIER.determineType(gene, event));
            }
        } else {
            String eventString = String.join(",", events);
            eventType.add(CLASSIFIER.determineType(gene, eventString));
        }
        return eventType;
    }
}
