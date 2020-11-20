package com.hartwig.hmftools.vicc.annotation;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class EventTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(EventTypeExtractor.class);

    @NotNull
    private static final EventClassifier CLASSIFIER = new EventClassifier(buildMatcherMap());

    public EventTypeExtractor() {
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

    @NotNull
    private static Map<EventType, EventMatcher> buildMatcherMap() {
        EventMatcher complexClassifier = new ComplexClassifier();
        EventMatcher combinedClassifier = new CombinedClassifier();
        EventMatcher fusionPairAndExonRangeClassifier = new FusionPairAndExonRangeClassifier();

        List<EventMatcher> firstTierEventMatchers =
                Lists.newArrayList(complexClassifier, combinedClassifier, fusionPairAndExonRangeClassifier);

        Map<EventType, EventMatcher> map = Maps.newHashMap();
        map.put(EventType.HOTSPOT, HotspotClassifier.create(firstTierEventMatchers));
        map.put(EventType.GENE_RANGE_CODON, GeneRangeCodonClassifier.create(firstTierEventMatchers));
        map.put(EventType.GENE_RANGE_EXON, GeneRangeExonClassifier.create(firstTierEventMatchers));
        map.put(EventType.FUSION_PAIR_AND_GENE_RANGE_EXON, fusionPairAndExonRangeClassifier);
        map.put(EventType.GENE_LEVEL, GeneLevelClassifier.create(firstTierEventMatchers));
        map.put(EventType.AMPLIFICATION, AmplificationClassifier.create(firstTierEventMatchers));
        map.put(EventType.DELETION, DeletionClassifier.create(firstTierEventMatchers));
        map.put(EventType.FUSION_PAIR, FusionPairClassifier.create(firstTierEventMatchers));
        map.put(EventType.PROMISCUOUS_FUSION, PromiscuousFusionClassifier.create(firstTierEventMatchers));
        map.put(EventType.SIGNATURE, SignatureClassifier.create(firstTierEventMatchers));
        map.put(EventType.COMBINED, combinedClassifier);
        map.put(EventType.COMPLEX, complexClassifier);

        return map;
    }
}
