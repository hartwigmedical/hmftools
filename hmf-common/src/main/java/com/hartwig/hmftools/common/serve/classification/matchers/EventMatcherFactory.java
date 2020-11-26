package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;
import com.hartwig.hmftools.common.serve.classification.MutationType;

import org.jetbrains.annotations.NotNull;

public final class EventMatcherFactory {

    private EventMatcherFactory() {
    }

    @NotNull
    public static Map<MutationType, EventMatcher> buildMatcherMap(@NotNull EventPreprocessor proteinAnnotationPreprocessor) {
        FusionPairMatcher fusionPairMatcher = new FusionPairMatcher();
        HotspotMatcher hotspotMatcher = new HotspotMatcher(proteinAnnotationPreprocessor, fusionPairMatcher);
        AmplificationMatcher amplificationMatcher = new AmplificationMatcher();

        EventMatcher complexClassifier = new ComplexMatcher();
        EventMatcher combinedClassifier = new CombinedMatcher(hotspotMatcher, fusionPairMatcher, amplificationMatcher);
        EventMatcher fusionPairAndExonRangeClassifier = new FusionPairAndExonRangeMatcher();

        List<EventMatcher> firstTierEventMatchers =
                Lists.newArrayList(complexClassifier, combinedClassifier, fusionPairAndExonRangeClassifier);

        Map<MutationType, EventMatcher> map = Maps.newHashMap();
        map.put(MutationType.HOTSPOT, withFirstTierMatchers(firstTierEventMatchers, hotspotMatcher));
        map.put(MutationType.GENE_RANGE_CODON,
                withFirstTierMatchers(firstTierEventMatchers, new GeneRangeCodonMatcher(proteinAnnotationPreprocessor)));
        map.put(MutationType.GENE_RANGE_EXON, withFirstTierMatchers(firstTierEventMatchers, new GeneRangeExonMatcher()));
        map.put(MutationType.FUSION_PAIR_AND_GENE_RANGE_EXON, fusionPairAndExonRangeClassifier);
        map.put(MutationType.GENE_LEVEL, withFirstTierMatchers(firstTierEventMatchers, new GeneLevelMatcher()));
        map.put(MutationType.AMPLIFICATION, withFirstTierMatchers(firstTierEventMatchers, amplificationMatcher));
        map.put(MutationType.DELETION, withFirstTierMatchers(firstTierEventMatchers, new DeletionMatcher()));
        map.put(MutationType.FUSION_PAIR, withFirstTierMatchers(firstTierEventMatchers, fusionPairMatcher));
        map.put(MutationType.PROMISCUOUS_FUSION,
                withFirstTierMatchers(firstTierEventMatchers, new PromiscuousFusionMatcher(fusionPairMatcher)));
        map.put(MutationType.SIGNATURE, withFirstTierMatchers(firstTierEventMatchers, new SignatureMatcher()));
        map.put(MutationType.COMBINED, combinedClassifier);
        map.put(MutationType.COMPLEX, complexClassifier);

        return map;
    }

    @NotNull
    private static EventMatcher withFirstTierMatchers(@NotNull List<EventMatcher> firstTierMatchers, @NotNull EventMatcher eventMatcher) {
        return new CompositeEventMatcher(firstTierMatchers, eventMatcher);
    }
}
