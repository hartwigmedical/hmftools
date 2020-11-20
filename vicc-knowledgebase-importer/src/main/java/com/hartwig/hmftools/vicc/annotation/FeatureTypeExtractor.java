package com.hartwig.hmftools.vicc.annotation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class FeatureTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(FeatureTypeExtractor.class);

    private FeatureTypeExtractor() {
    }

    @NotNull
    public static FeatureType extractType(@NotNull Feature feature) {
        String gene = feature.geneSymbol();
        if (gene == null) {
            LOGGER.debug("Skipping extraction for '{}' since gene is missing", feature.name());
            return FeatureType.UNKNOWN;
        } else {
            return extractType(gene, feature.name());
        }
    }

    @NotNull
    private static FeatureType extractType(@NotNull String gene, @NotNull String event) {
        Map<FeatureType, Boolean> evaluations = Maps.newHashMap();

        evaluations.put(FeatureType.HOTSPOT, HotspotClassifier.isHotspot(event));
        evaluations.put(FeatureType.GENE_RANGE_CODON, GeneRangeClassifier.isGeneRangeCodonEvent(event));
        evaluations.put(FeatureType.GENE_RANGE_EXON, GeneRangeClassifier.isGeneRangeExonEvent(gene, event));
        evaluations.put(FeatureType.GENE_LEVEL, GeneRangeClassifier.isGeneLevelEvent(gene, event));
        evaluations.put(FeatureType.AMPLIFICATION, CopyNumberClassifier.isAmplification(gene, event));
        evaluations.put(FeatureType.DELETION, CopyNumberClassifier.isDeletion(gene, event));
        evaluations.put(FeatureType.FUSION_PAIR, FusionClassifier.isFusionPair(gene, event));
        evaluations.put(FeatureType.PROMISCUOUS_FUSION, FusionClassifier.isPromiscuousFusion(gene, event));
        evaluations.put(FeatureType.FUSION_PAIR_AND_GENE_RANGE_EXON, CombinedClassifier.isFusionPairAndGeneRangeExon(gene, event));
        evaluations.put(FeatureType.SIGNATURE, SignatureClassifier.isSignature(event));
        evaluations.put(FeatureType.COMBINED, CombinedClassifier.isCombinedEvent(gene, event));
        evaluations.put(FeatureType.COMPLEX, ComplexClassifier.isComplexEvent(gene, event));

        Set<FeatureType> positiveTypes = Sets.newHashSet();
        for (Map.Entry<FeatureType, Boolean> evaluation : evaluations.entrySet()) {
            if (evaluation.getValue()) {
                positiveTypes.add(evaluation.getKey());
            }
        }

        if (positiveTypes.size() > 1) {
            LOGGER.warn("More than one type evaluated to true for '{}' on '{}': {}", event, gene, positiveTypes);
        } else if (positiveTypes.size() == 1) {
            return positiveTypes.iterator().next();
        }

        return FeatureType.UNKNOWN;
    }
}
