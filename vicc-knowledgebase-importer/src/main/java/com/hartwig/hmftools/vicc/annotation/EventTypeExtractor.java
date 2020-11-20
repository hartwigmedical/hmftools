package com.hartwig.hmftools.vicc.annotation;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class EventTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(EventTypeExtractor.class);

    @NotNull
    private static final Map<EventType, EventClassifier> MATCHERS = buildMatcherMap();

    public EventTypeExtractor() {
    }

    @NotNull
    public static EventType extractType(@NotNull Feature feature) {
        String gene = feature.geneSymbol();
        if (gene == null) {
            LOGGER.debug("Skipping extraction for '{}' since gene is missing", feature.name());
            return EventType.UNKNOWN;
        } else {
            return extractType(gene, feature.name());
        }
    }

    @NotNull
    private static EventType extractType(@NotNull String gene, @NotNull String event) {
        Map<EventType, Boolean> evaluations = Maps.newHashMap();

        for (Map.Entry<EventType, EventClassifier> entry : MATCHERS.entrySet()){
            evaluations.put(entry.getKey(), entry.getValue().matches(gene, event));
        }

        Set<EventType> positiveTypes = Sets.newHashSet();
        for (Map.Entry<EventType, Boolean> evaluation : evaluations.entrySet()) {
            if (evaluation.getValue()) {
                positiveTypes.add(evaluation.getKey());
            }
        }

        if (positiveTypes.size() > 1) {
            LOGGER.warn("More than one type evaluated to true for '{}' on '{}': {}", event, gene, positiveTypes);
        } else if (positiveTypes.size() == 1) {
            return positiveTypes.iterator().next();
        }

        return EventType.UNKNOWN;
    }

    @NotNull
    private static Map<EventType, EventClassifier> buildMatcherMap() {
        EventClassifier complexClassifier = new ComplexClassifier();
        EventClassifier combinedClassifier = new CombinedClassifier();
        EventClassifier fusionPairAndExonRangeClassifier = new FusionPairAndExonRangeClassifier();

        List<EventClassifier> firstTierEventClassifiers =
                Lists.newArrayList(complexClassifier, combinedClassifier, fusionPairAndExonRangeClassifier);

        Map<EventType, EventClassifier> map = Maps.newHashMap();
        map.put(EventType.HOTSPOT, HotspotClassifier.create(firstTierEventClassifiers));
        map.put(EventType.GENE_RANGE_CODON, GeneRangeCodonClassifier.create(firstTierEventClassifiers));
        map.put(EventType.GENE_RANGE_EXON, GeneRangeExonClassifier.create(firstTierEventClassifiers));
        map.put(EventType.FUSION_PAIR_AND_GENE_RANGE_EXON, fusionPairAndExonRangeClassifier);
        map.put(EventType.GENE_LEVEL, GeneLevelClassifier.create(firstTierEventClassifiers));
        map.put(EventType.AMPLIFICATION, AmplificationClassifier.create(firstTierEventClassifiers));
        map.put(EventType.DELETION, DeletionClassifier.create(firstTierEventClassifiers));
        map.put(EventType.FUSION_PAIR, FusionPairClassifier.create(firstTierEventClassifiers));
        map.put(EventType.PROMISCUOUS_FUSION, PromiscuousFusionClassifier.create(firstTierEventClassifiers));
        map.put(EventType.SIGNATURE, SignatureClassifier.create(firstTierEventClassifiers));
        map.put(EventType.COMBINED, combinedClassifier);
        map.put(EventType.COMPLEX, complexClassifier);

        return map;
    }
}
