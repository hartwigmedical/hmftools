package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.common.serve.classification.MutationType;

import org.jetbrains.annotations.NotNull;

public final class EventMatcherFactory {

    private EventMatcherFactory() {
    }

    @NotNull
    public static Map<MutationType, EventMatcher> buildMatcherMap(@NotNull EventClassifierConfig config) {
        FusionPairMatcher fusionPairMatcher = new FusionPairMatcher(config.exonicDelDupFusionEvents(), config.fusionPairEventsToSkip());
        HotspotMatcher hotspotMatcher = new HotspotMatcher(config.proteinAnnotationExtractor(), fusionPairMatcher);
        AmplificationMatcher amplificationMatcher =
                new AmplificationMatcher(config.amplificationKeywords(), config.amplificationKeyPhrases());

        EventMatcher complexClassifier = new ComplexMatcher(config.complexEventsPerGene());
        EventMatcher combinedClassifier =
                new CombinedMatcher(config.combinedEventsPerGene(), hotspotMatcher, fusionPairMatcher, amplificationMatcher);
        EventMatcher fusionPairAndExonRangeClassifier = new FusionPairAndExonRangeMatcher(config.fusionPairAndExonRangesPerGene());

        List<EventMatcher> firstTierEventMatchers =
                Lists.newArrayList(complexClassifier, combinedClassifier, fusionPairAndExonRangeClassifier);

        EventMatcher geneRangeCodonMatcher = new GeneRangeCodonMatcher(config.proteinAnnotationExtractor());
        EventMatcher geneRangeExonMatcher =
                new GeneRangeExonMatcher(config.exonKeyword(), config.exonRangeEvents(), config.exonRangeKeywords());
        EventMatcher geneLevelMatcher = new GeneLevelMatcher(config.exonKeyword(),
                config.genericGeneLevelKeywords(),
                config.activatingGeneLevelKeywords(),
                config.inactivatingGeneLevelKeywords());
        EventMatcher deletionMatcher =
                new DeletionMatcher(config.deletionKeywords(), config.deletionKeyPhrases(), config.deletionKeywordsToSkip());
        EventMatcher promiscuousFusionMatcher = new PromiscuousFusionMatcher(config.promiscuousFusionKeywords(), fusionPairMatcher);
        EventMatcher signatureMatcher = new SignatureMatcher(config.signatureEvents());

        Map<MutationType, EventMatcher> map = Maps.newHashMap();
        map.put(MutationType.HOTSPOT, withFirstTierMatchers(firstTierEventMatchers, hotspotMatcher));
        map.put(MutationType.GENE_RANGE_CODON, withFirstTierMatchers(firstTierEventMatchers, geneRangeCodonMatcher));
        map.put(MutationType.GENE_RANGE_EXON, withFirstTierMatchers(firstTierEventMatchers, geneRangeExonMatcher));
        map.put(MutationType.FUSION_PAIR_AND_GENE_RANGE_EXON, fusionPairAndExonRangeClassifier);
        map.put(MutationType.GENE_LEVEL, withFirstTierMatchers(firstTierEventMatchers, geneLevelMatcher));
        map.put(MutationType.AMPLIFICATION, withFirstTierMatchers(firstTierEventMatchers, amplificationMatcher));
        map.put(MutationType.DELETION, withFirstTierMatchers(firstTierEventMatchers, deletionMatcher));
        map.put(MutationType.FUSION_PAIR, withFirstTierMatchers(firstTierEventMatchers, fusionPairMatcher));
        map.put(MutationType.PROMISCUOUS_FUSION, withFirstTierMatchers(firstTierEventMatchers, promiscuousFusionMatcher));
        map.put(MutationType.SIGNATURE, withFirstTierMatchers(firstTierEventMatchers, signatureMatcher));
        map.put(MutationType.COMBINED, combinedClassifier);
        map.put(MutationType.COMPLEX, complexClassifier);

        return map;
    }

    @NotNull
    private static EventMatcher withFirstTierMatchers(@NotNull List<EventMatcher> firstTierMatchers, @NotNull EventMatcher eventMatcher) {
        return new CompositeEventMatcher(firstTierMatchers, eventMatcher);
    }
}
