package com.hartwig.hmftools.vicc.annotation;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class MutationTypeExtractor {

    private static final Logger LOGGER = LogManager.getLogger(MutationTypeExtractor.class);

    @NotNull
    private static final EventClassifier CLASSIFIER = new EventClassifier(buildMatcherMap());

    public MutationTypeExtractor() {
    }

    @NotNull
    public static MutationType extractType(@NotNull Feature feature) {
        String gene = feature.geneSymbol();
        if (gene == null) {
            LOGGER.debug("Skipping extraction for '{}' since gene is missing", feature.name());
            return MutationType.UNKNOWN;
        } else {
            return CLASSIFIER.determineType(gene, feature.name());
        }
    }

    @NotNull
    private static Map<MutationType, EventMatcher> buildMatcherMap() {
        EventMatcher complexClassifier = new ComplexClassifier();
        EventMatcher combinedClassifier = new CombinedClassifier();
        EventMatcher fusionPairAndExonRangeClassifier = new FusionPairAndExonRangeClassifier();

        List<EventMatcher> firstTierEventMatchers =
                Lists.newArrayList(complexClassifier, combinedClassifier, fusionPairAndExonRangeClassifier);

        Map<MutationType, EventMatcher> map = Maps.newHashMap();
        map.put(MutationType.HOTSPOT, HotspotClassifier.create(firstTierEventMatchers));
        map.put(MutationType.GENE_RANGE_CODON, GeneRangeCodonClassifier.create(firstTierEventMatchers));
        map.put(MutationType.GENE_RANGE_EXON, GeneRangeExonClassifier.create(firstTierEventMatchers));
        map.put(MutationType.FUSION_PAIR_AND_GENE_RANGE_EXON, fusionPairAndExonRangeClassifier);
        map.put(MutationType.GENE_LEVEL, GeneLevelClassifier.create(firstTierEventMatchers));
        map.put(MutationType.AMPLIFICATION, AmplificationClassifier.create(firstTierEventMatchers));
        map.put(MutationType.DELETION, DeletionClassifier.create(firstTierEventMatchers));
        map.put(MutationType.FUSION_PAIR, FusionPairClassifier.create(firstTierEventMatchers));
        map.put(MutationType.PROMISCUOUS_FUSION, PromiscuousFusionClassifier.create(firstTierEventMatchers));
        map.put(MutationType.SIGNATURE, SignatureClassifier.create(firstTierEventMatchers));
        map.put(MutationType.COMBINED, combinedClassifier);
        map.put(MutationType.COMPLEX, complexClassifier);

        return map;
    }
}
