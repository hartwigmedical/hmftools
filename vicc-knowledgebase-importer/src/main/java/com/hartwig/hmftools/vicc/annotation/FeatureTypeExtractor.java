package com.hartwig.hmftools.vicc.annotation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FeatureTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(FeatureTypeExtractor.class);

    private FeatureTypeExtractor() {
    }

    @NotNull
    public static FeatureType extractType(@NotNull Feature feature) {
        return extractType(feature.name(), feature.geneSymbol(), feature.biomarkerType(), feature.provenanceRule());
    }

    @NotNull
    public static FeatureType extractType(@NotNull String featureName, @Nullable String gene, @Nullable String biomarkerType,
            @Nullable String provenanceRule) {
        Map<FeatureType, Boolean> evaluations = Maps.newHashMap();

        evaluations.put(FeatureType.HOTSPOT, HotspotClassifier.isHotspot(featureName));
        evaluations.put(FeatureType.GENE_RANGE_CODON, GeneRangeClassifier.isGeneRangeCodonEvent(featureName));
        evaluations.put(FeatureType.GENE_RANGE_EXON, GeneRangeClassifier.isGeneRangeExonEvent(featureName, gene));
        evaluations.put(FeatureType.GENE_LEVEL, GeneRangeClassifier.isGeneLevelEvent(featureName, provenanceRule));
        evaluations.put(FeatureType.AMPLIFICATION, CopyNumberClassifier.isAmplification(featureName, biomarkerType));
        evaluations.put(FeatureType.DELETION, CopyNumberClassifier.isDeletion(featureName, biomarkerType));
        evaluations.put(FeatureType.FUSION_PAIR, FusionClassifier.isFusionPair(featureName, gene, biomarkerType));
        evaluations.put(FeatureType.FUSION_PROMISCUOUS, FusionClassifier.isPromiscuousFusion(featureName, gene, biomarkerType));
        evaluations.put(FeatureType.FUSION_PAIR_AND_GENE_RANGE_EXON, CombinedClassifier.isFusionPairAndGeneRangeExon(featureName, gene));
        evaluations.put(FeatureType.SIGNATURE, SignatureClassifier.isSignature(featureName));
        evaluations.put(FeatureType.COMBINED, CombinedClassifier.isCombinedEvent(featureName));

        Set<FeatureType> positiveTypes = Sets.newHashSet();
        for (Map.Entry<FeatureType, Boolean> evaluation : evaluations.entrySet()) {
            if (evaluation.getValue()) {
                positiveTypes.add(evaluation.getKey());
            }
        }

        if (positiveTypes.size() > 1) {
            LOGGER.warn("More than one type evaluated to true for '{}': {}", featureName, positiveTypes);
        } else if (positiveTypes.size() == 1) {
            return positiveTypes.iterator().next();
        }

        return FeatureType.UNKNOWN;
    }
}
