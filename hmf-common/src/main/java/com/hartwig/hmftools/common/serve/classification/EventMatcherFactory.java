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
        EventMatcher complexClassifier = new ComplexMatcher();
        EventMatcher combinedClassifier = new CombinedMatcher();
        EventMatcher fusionPairAndExonRangeClassifier = new FusionPairAndExonRangeMatcher();

        List<EventMatcher> firstTierEventMatchers =
                Lists.newArrayList(complexClassifier, combinedClassifier, fusionPairAndExonRangeClassifier);

        Map<MutationType, EventMatcher> map = Maps.newHashMap();
        map.put(MutationType.HOTSPOT, HotspotMatcher.create(firstTierEventMatchers, proteinAnnotationPreprocessor));
        map.put(MutationType.GENE_RANGE_CODON, GeneRangeCodonMatcher.create(firstTierEventMatchers, proteinAnnotationPreprocessor));
        map.put(MutationType.GENE_RANGE_EXON, GeneRangeExonMatcher.create(firstTierEventMatchers));
        map.put(MutationType.FUSION_PAIR_AND_GENE_RANGE_EXON, fusionPairAndExonRangeClassifier);
        map.put(MutationType.GENE_LEVEL, GeneLevelMatcher.create(firstTierEventMatchers));
        map.put(MutationType.AMPLIFICATION, AmplificationMatcher.create(firstTierEventMatchers));
        map.put(MutationType.DELETION, DeletionMatcher.create(firstTierEventMatchers));
        map.put(MutationType.FUSION_PAIR, FusionPairMatcher.create(firstTierEventMatchers));
        map.put(MutationType.PROMISCUOUS_FUSION, PromiscuousFusionMatcher.create(firstTierEventMatchers));
        map.put(MutationType.SIGNATURE, SignatureMatcher.create(firstTierEventMatchers));
        map.put(MutationType.COMBINED, combinedClassifier);
        map.put(MutationType.COMPLEX, complexClassifier);

        return map;
    }
}
