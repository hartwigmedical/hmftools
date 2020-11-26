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

        ComplexMatcher complexMatcher = new ComplexMatcher(config.complexEventsPerGene());
        CombinedMatcher combinedMatcher =
                new CombinedMatcher(config.combinedEventsPerGene(), hotspotMatcher, fusionPairMatcher, amplificationMatcher);
        FusionPairAndExonMatcher fusionPairAndExonMatcher = new FusionPairAndExonMatcher(config.fusionPairAndExonsPerGene());

        List<EventMatcher> firstTierEventMatchers = Lists.newArrayList(complexMatcher, combinedMatcher, fusionPairAndExonMatcher);

        CodonMatcher codonMatcher = new CodonMatcher(config.proteinAnnotationExtractor());
        ExonMatcher exonMatcher = new ExonMatcher(config.exonIdentifiers(), config.exonKeywords(), config.specificExonEvents());
        GeneLevelMatcher geneLevelMatcher = new GeneLevelMatcher(config.geneLevelBlacklistKeyPhrases(),
                config.genericGeneLevelKeyPhrases(),
                config.activatingGeneLevelKeyPhrases(),
                config.inactivatingGeneLevelKeyPhrases());
        DeletionMatcher deletionMatcher =
                new DeletionMatcher(config.deletionKeywords(), config.deletionKeyPhrases(), config.deletionKeyPhrasesToSkip());
        PromiscuousFusionMatcher promiscuousFusionMatcher =
                new PromiscuousFusionMatcher(config.promiscuousFusionKeyPhrases(), fusionPairMatcher);
        SignatureMatcher signatureMatcher = new SignatureMatcher(config.signatureEvents());

        Map<MutationType, EventMatcher> map = Maps.newHashMap();
        map.put(MutationType.HOTSPOT, withFirstTierMatchers(firstTierEventMatchers, hotspotMatcher));
        map.put(MutationType.CODON, withFirstTierMatchers(firstTierEventMatchers, codonMatcher));
        map.put(MutationType.EXON, withFirstTierMatchers(firstTierEventMatchers, exonMatcher));
        map.put(MutationType.FUSION_PAIR_AND_EXON, fusionPairAndExonMatcher);
        map.put(MutationType.GENE_LEVEL, withFirstTierMatchers(firstTierEventMatchers, geneLevelMatcher));
        map.put(MutationType.AMPLIFICATION, withFirstTierMatchers(firstTierEventMatchers, amplificationMatcher));
        map.put(MutationType.DELETION, withFirstTierMatchers(firstTierEventMatchers, deletionMatcher));
        map.put(MutationType.FUSION_PAIR, withFirstTierMatchers(firstTierEventMatchers, fusionPairMatcher));
        map.put(MutationType.PROMISCUOUS_FUSION, withFirstTierMatchers(firstTierEventMatchers, promiscuousFusionMatcher));
        map.put(MutationType.SIGNATURE, withFirstTierMatchers(firstTierEventMatchers, signatureMatcher));
        map.put(MutationType.COMBINED, combinedMatcher);
        map.put(MutationType.COMPLEX, complexMatcher);

        return map;
    }

    @NotNull
    private static EventMatcher withFirstTierMatchers(@NotNull List<EventMatcher> firstTierMatchers, @NotNull EventMatcher eventMatcher) {
        return new CompositeEventMatcher(firstTierMatchers, eventMatcher);
    }
}
