package com.hartwig.hmftools.vicc.annotation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class EventTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(EventTypeExtractor.class);

    private EventTypeExtractor() {
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

        evaluations.put(EventType.HOTSPOT, HotspotClassifier.isHotspot(event));
        evaluations.put(EventType.GENE_RANGE_CODON, GeneRangeClassifier.isGeneRangeCodonEvent(event));
        evaluations.put(EventType.GENE_RANGE_EXON, GeneRangeClassifier.isGeneRangeExonEvent(gene, event));
        evaluations.put(EventType.GENE_LEVEL, GeneRangeClassifier.isGeneLevelEvent(gene, event));
        evaluations.put(EventType.AMPLIFICATION, CopyNumberClassifier.isAmplification(gene, event));
        evaluations.put(EventType.DELETION, CopyNumberClassifier.isDeletion(gene, event));
        evaluations.put(EventType.FUSION_PAIR, FusionClassifier.isFusionPair(gene, event));
        evaluations.put(EventType.PROMISCUOUS_FUSION, FusionClassifier.isPromiscuousFusion(gene, event));
        evaluations.put(EventType.FUSION_PAIR_AND_GENE_RANGE_EXON, CombinedClassifier.isFusionPairAndGeneRangeExon(gene, event));
        evaluations.put(EventType.SIGNATURE, SignatureClassifier.isSignature(event));
        evaluations.put(EventType.COMBINED, CombinedClassifier.isCombinedEvent(gene, event));
        evaluations.put(EventType.COMPLEX, ComplexClassifier.isComplexEvent(gene, event));

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
}
