package com.hartwig.hmftools.common.serve.classification;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class EventMatcherFactory {

    private EventMatcherFactory() {
    }

    @NotNull
    public static Map<MutationType, EventMatcher> buildMatcherMap(@NotNull EventPreprocessor proteinAnnotationPreprocessor) {
        EventMatcher complexClassifier = new ComplexClassifier();
        EventMatcher combinedClassifier = new CombinedClassifier();
        EventMatcher fusionPairAndExonRangeClassifier = new FusionPairAndExonRangeClassifier();

        List<EventMatcher> firstTierEventMatchers =
                Lists.newArrayList(complexClassifier, combinedClassifier, fusionPairAndExonRangeClassifier);

        Map<MutationType, EventMatcher> map = Maps.newHashMap();
        map.put(MutationType.HOTSPOT, HotspotClassifier.create(firstTierEventMatchers, proteinAnnotationPreprocessor));
        map.put(MutationType.GENE_RANGE_CODON, GeneRangeCodonClassifier.create(firstTierEventMatchers, proteinAnnotationPreprocessor));
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
