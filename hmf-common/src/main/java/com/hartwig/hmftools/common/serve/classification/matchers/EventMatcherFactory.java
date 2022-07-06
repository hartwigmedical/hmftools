package com.hartwig.hmftools.common.serve.classification.matchers;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.jetbrains.annotations.NotNull;

public final class EventMatcherFactory {

    private EventMatcherFactory() {
    }

    @NotNull
    public static Map<EventType, EventMatcher> buildMatcherMap(@NotNull EventClassifierConfig config) {
        FusionPairMatcher fusionPairMatcher = new FusionPairMatcher(config.exonicDelDupFusionKeyPhrases(),
                config.exonicDelDupFusionEvents(),
                config.fusionPairEventsToSkip());
        PromiscuousFusionMatcher promiscuousFusionMatcher =
                new PromiscuousFusionMatcher(config.promiscuousFusionKeyPhrases(), fusionPairMatcher);

        HotspotMatcher hotspotMatcher = new HotspotMatcher(config.proteinAnnotationExtractor(), fusionPairMatcher);
        CodonMatcher codonMatcher = new CodonMatcher(config.proteinAnnotationExtractor());
        ExonMatcher exonMatcher = new ExonMatcher(config.exonIdentifiers(),
                config.exonKeywords(),
                config.exonBlacklistKeyPhrases(),
                config.specificExonEvents());
        GeneLevelMatcher geneLevelMatcher = new GeneLevelMatcher(config.geneLevelBlacklistKeyPhrases(),
                config.genericGeneLevelKeyPhrases(),
                config.activatingGeneLevelKeyPhrases(),
                config.inactivatingGeneLevelKeyPhrases());

        WildTypeMatcher wildTypeMatcher = new WildTypeMatcher(config.wildTypeKeyPhrases());

        FusionPairAndExonMatcher fusionPairAndExonMatcher = new FusionPairAndExonMatcher(config.fusionPairAndExonsPerGene());
        AmplificationMatcher amplificationMatcher = new AmplificationMatcher(config.amplificationKeywords(),
                config.amplificationKeyPhrases());
        DeletionMatcher deletionMatcher = new DeletionMatcher(config.deletionBlacklistKeyPhrases(),
                config.deletionKeywords(),
                config.deletionKeyPhrases());

        OverExpressionMatcher overExpressionMatcher = new OverExpressionMatcher(config.overexpressionKeywords(),
                config.overexpressionKeyPhrases());
        UnderExpressionMatcher underExpressionMatcher = new UnderExpressionMatcher(
                config.underexpressionKeywords(),
                config.underexpressionKeyPhrases());

        CharacteristicMatcher characteristicMatcher = new CharacteristicMatcher(allCharacteristicKeyPhrases(config));
        HlaMatcher hlaMatcher = new HlaMatcher(allHlaKeyPhrases(config));

        ComplexMatcher complexMatcher = new ComplexMatcher(hotspotMatcher, config.complexEventsPerGene());
        CombinedMatcher combinedMatcher =
                new CombinedMatcher(config.combinedEventsPerGene(), hotspotMatcher, fusionPairMatcher, amplificationMatcher);

        List<EventMatcher> firstTierEventMatchers = Lists.newArrayList(complexMatcher, combinedMatcher, fusionPairAndExonMatcher);

        Map<EventType, EventMatcher> map = Maps.newHashMap();
        map.put(EventType.HOTSPOT, withFirstTierMatchers(firstTierEventMatchers, hotspotMatcher));
        map.put(EventType.CODON, withFirstTierMatchers(firstTierEventMatchers, codonMatcher));
        map.put(EventType.EXON, withFirstTierMatchers(firstTierEventMatchers, exonMatcher));
        map.put(EventType.FUSION_PAIR_AND_EXON, fusionPairAndExonMatcher);
        map.put(EventType.GENE_LEVEL, withFirstTierMatchers(firstTierEventMatchers, geneLevelMatcher));
        map.put(EventType.WILD_TYPE, withFirstTierMatchers(firstTierEventMatchers, wildTypeMatcher));
        map.put(EventType.AMPLIFICATION, withFirstTierMatchers(firstTierEventMatchers, amplificationMatcher));
        map.put(EventType.OVER_EXPRESSION, withFirstTierMatchers(firstTierEventMatchers, overExpressionMatcher));
        map.put(EventType.DELETION, withFirstTierMatchers(firstTierEventMatchers, deletionMatcher));
        map.put(EventType.UNDER_EXPRESSION, withFirstTierMatchers(firstTierEventMatchers, underExpressionMatcher));
        map.put(EventType.FUSION_PAIR, withFirstTierMatchers(firstTierEventMatchers, fusionPairMatcher));
        map.put(EventType.PROMISCUOUS_FUSION, withFirstTierMatchers(firstTierEventMatchers, promiscuousFusionMatcher));
        map.put(EventType.CHARACTERISTIC, withFirstTierMatchers(firstTierEventMatchers, characteristicMatcher));
        map.put(EventType.IMMUNO_HLA, withFirstTierMatchers(firstTierEventMatchers, hlaMatcher));
        map.put(EventType.COMBINED, combinedMatcher);
        map.put(EventType.COMPLEX, complexMatcher);

        return map;
    }

    @NotNull
    private static Set<String> allCharacteristicKeyPhrases(@NotNull EventClassifierConfig config) {
        Set<String> tumorCharacteristics = Sets.newHashSet();
        tumorCharacteristics.addAll(config.microsatelliteUnstableKeyPhrases());
        tumorCharacteristics.addAll(config.microsatelliteStableKeyPhrases());
        tumorCharacteristics.addAll(config.highTumorMutationalLoadKeyPhrases());
        tumorCharacteristics.addAll(config.lowTumorMutationalLoadKeyPhrases());
        tumorCharacteristics.addAll(config.highTumorMutationalBurdenKeyPhrases());
        tumorCharacteristics.addAll(config.lowTumorMutationalBurdenKeyPhrases());
        tumorCharacteristics.addAll(config.hrDeficiencyKeyPhrases());
        tumorCharacteristics.addAll(config.hpvPositiveEvents());
        tumorCharacteristics.addAll(config.ebvPositiveEvents());
        return tumorCharacteristics;
    }

    @NotNull
    private static Set<String> allHlaKeyPhrases(@NotNull EventClassifierConfig config) {
        Set<String> hla = Sets.newHashSet();
        hla.addAll(config.hlaKeyPhrases());
        return hla;
    }

    @NotNull
    private static EventMatcher withFirstTierMatchers(@NotNull List<EventMatcher> firstTierMatchers, @NotNull EventMatcher eventMatcher) {
        return new CompositeEventMatcher(firstTierMatchers, eventMatcher);
    }
}