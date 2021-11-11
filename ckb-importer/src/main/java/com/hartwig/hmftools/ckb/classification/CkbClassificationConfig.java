package com.hartwig.hmftools.ckb.classification;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.common.serve.classification.ImmutableEventClassifierConfig;

import org.jetbrains.annotations.NotNull;

public class CkbClassificationConfig {

    private static final Set<String> EXON_IDENTIFIERS = exonIdentifiers();
    private static final Set<String> EXON_KEYWORDS = exonKeywords();
    private static final Set<String> EXON_BLACKLIST_KEY_PHRASES = exonBlacklistKeyPhrases();
    private static final Set<String> SPECIFIC_EXON_EVENTS = specificExonEvents();
    private static final Map<String, Set<String>> FUSION_PAIR_AND_EXONS_PER_GENE = fusionPairAndExonsPerGene();
    private static final Set<String> GENE_LEVEL_BLACKLIST_KEY_PHRASES = geneLevelBlacklistKeyPhrases();
    private static final Set<String> GENERIC_GENE_LEVEL_KEY_PHRASES = genericGeneLevelKeyPhrases();
    private static final Set<String> ACTIVATING_GENE_LEVEL_KEY_PHRASES = activatingGeneLevelKeyPhrases();
    private static final Set<String> INACTIVATING_GENE_LEVEL_KEY_PHRASES = inactivatingGeneLevelKeyPhrases();
    private static final Set<String> AMPLIFICATION_KEYWORDS = amplificationKeywords();
    private static final Set<String> AMPLIFICATION_KEY_PHRASES = amplificationKeyPhrases();
    private static final Set<String> DELETION_BLACKLIST_KEY_PHRASES = deletionBlacklistKeyPhrases();
    private static final Set<String> DELETION_KEYWORDS = deletionKeywords();
    private static final Set<String> DELETION_KEY_PHRASES = deletionKeyPhrases();
    private static final Set<String> EXONIC_DEL_DUP_FUSION_KEY_PHRASES = exonicDelDupFusionKeyPhrases();
    private static final Set<String> EXONIC_DEL_DUP_FUSION_EVENTS = exonicDelDupFusionEvents();
    private static final Set<String> FUSION_PAIR_EVENTS_TO_SKIP = fusionPairEventsToSkip();
    private static final Set<String> PROMISCUOUS_FUSION_KEY_PHRASES = promiscuousFusionKeyPhrases();
    private static final Set<String> MICROSATELLITE_UNSTABLE_EVENTS = microsatelliteUnstableEvents();
    private static final Set<String> MICROSATELLITE_STABLE_EVENTS = microsatelliteStableEvents();
    private static final Set<String> HIGH_TUMOR_MUTATIONAL_LOAD_EVENTS = highTumorMutationalLoadEvents();
    private static final Set<String> LOW_TUMOR_MUTATIONAL_LOAD_EVENTS = lowTumorMutationalLoadEvents();
    private static final Set<String> HR_DEFICIENCY_EVENTS = hrDeficiencyEvents();
    private static final Set<String> HLA_EVENTS = hlaEvents();
    private static final Set<String> HPV_POSITIVE_EVENTS = hpvPositiveEvents();
    private static final Set<String> EBV_POSITIVE_EVENTS = ebvPositiveEvents();
    private static final Map<String, Set<String>> COMBINED_EVENTS_PER_GENE = combinedEventsPerGene();
    private static final Map<String, Set<String>> COMPLEX_EVENTS_PER_GENE = complexEventsPerGene();

    private CkbClassificationConfig() {
    }

    @NotNull
    public static EventClassifierConfig build() {
        return ImmutableEventClassifierConfig.builder()
                .proteinAnnotationExtractor(new CkbProteinAnnotationExtractor())
                .exonIdentifiers(EXON_IDENTIFIERS)
                .exonKeywords(EXON_KEYWORDS)
                .exonBlacklistKeyPhrases(EXON_BLACKLIST_KEY_PHRASES)
                .specificExonEvents(SPECIFIC_EXON_EVENTS)
                .fusionPairAndExonsPerGene(FUSION_PAIR_AND_EXONS_PER_GENE)
                .geneLevelBlacklistKeyPhrases(GENE_LEVEL_BLACKLIST_KEY_PHRASES)
                .genericGeneLevelKeyPhrases(GENERIC_GENE_LEVEL_KEY_PHRASES)
                .activatingGeneLevelKeyPhrases(ACTIVATING_GENE_LEVEL_KEY_PHRASES)
                .inactivatingGeneLevelKeyPhrases(INACTIVATING_GENE_LEVEL_KEY_PHRASES)
                .amplificationKeywords(AMPLIFICATION_KEYWORDS)
                .amplificationKeyPhrases(AMPLIFICATION_KEY_PHRASES)
                .deletionBlacklistKeyPhrases(DELETION_BLACKLIST_KEY_PHRASES)
                .deletionKeywords(DELETION_KEYWORDS)
                .deletionKeyPhrases(DELETION_KEY_PHRASES)
                .exonicDelDupFusionKeyPhrases(EXONIC_DEL_DUP_FUSION_KEY_PHRASES)
                .exonicDelDupFusionEvents(EXONIC_DEL_DUP_FUSION_EVENTS)
                .fusionPairEventsToSkip(FUSION_PAIR_EVENTS_TO_SKIP)
                .promiscuousFusionKeyPhrases(PROMISCUOUS_FUSION_KEY_PHRASES)
                .microsatelliteUnstableEvents(MICROSATELLITE_UNSTABLE_EVENTS)
                .microsatelliteStableEvents(MICROSATELLITE_STABLE_EVENTS)
                .highTumorMutationalLoadEvents(HIGH_TUMOR_MUTATIONAL_LOAD_EVENTS)
                .lowTumorMutationalLoadEvents(LOW_TUMOR_MUTATIONAL_LOAD_EVENTS)
                .hrDeficiencyEvents(HR_DEFICIENCY_EVENTS)
                .hlaEvents(HLA_EVENTS)
                .hpvPositiveEvents(HPV_POSITIVE_EVENTS)
                .ebvPositiveEvents(EBV_POSITIVE_EVENTS)
                .combinedEventsPerGene(COMBINED_EVENTS_PER_GENE)
                .complexEventsPerGene(COMPLEX_EVENTS_PER_GENE)
                .build();
    }

    @NotNull
    private static Set<String> exonIdentifiers() {
        Set<String> set = Sets.newHashSet();
        set.add("exon");
        return set;
    }

    @NotNull
    private static Set<String> exonKeywords() {
        Set<String> set = Sets.newHashSet();
        set.add("exon");
        return set;
    }

    @NotNull
    private static Set<String> exonBlacklistKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("del exon");
        return set;
    }

    @NotNull
    private static Set<String> specificExonEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Map<String, Set<String>> fusionPairAndExonsPerGene() {
        Map<String, Set<String>> map = Maps.newHashMap();
        map.put("KIT", Sets.newHashSet("exon 11 del", "exon 11"));

        return map;
    }

    @NotNull
    private static Set<String> geneLevelBlacklistKeyPhrases() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> genericGeneLevelKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("mutant");
        return set;
    }

    @NotNull
    private static Set<String> activatingGeneLevelKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("act mut");
        return set;
    }

    @NotNull
    private static Set<String> inactivatingGeneLevelKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("inact mut");
        set.add("negative");
        return set;
    }

    @NotNull
    private static Set<String> amplificationKeywords() {
        Set<String> set = Sets.newHashSet();
        set.add("amp");
        return set;
    }

    @NotNull
    private static Set<String> amplificationKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("over exp");
        return set;
    }

    @NotNull
    private static Set<String> deletionBlacklistKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("exon");
        return set;
    }

    @NotNull
    private static Set<String> deletionKeywords() {
        Set<String> set = Sets.newHashSet();
        set.add("loss");
        set.add("del");
        return set;
    }

    @NotNull
    private static Set<String> deletionKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("dec exp");
        return set;
    }

    @NotNull
    private static Set<String> exonicDelDupFusionKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("del exon");
        return set;
    }

    @NotNull
    private static Set<String> exonicDelDupFusionEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> fusionPairEventsToSkip() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> promiscuousFusionKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("fusion promiscuous");
        set.add("rearrange");
        return set;
    }

    @NotNull
    private static Set<String> microsatelliteUnstableEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("MSI high");
        return set;
    }

    @NotNull
    private static Set<String> microsatelliteStableEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("MSI low");
        set.add("MSI neg");
        return set;
    }

    @NotNull
    private static Set<String> highTumorMutationalLoadEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("TMB high");
        return set;
    }

    @NotNull
    private static Set<String> lowTumorMutationalLoadEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("TMB low");
        return set;
    }

    @NotNull
    private static Set<String> hrDeficiencyEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> hlaEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> hpvPositiveEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> ebvPositiveEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Map<String, Set<String>> combinedEventsPerGene() {
        return Maps.newHashMap();
    }

    @NotNull
    private static Map<String, Set<String>> complexEventsPerGene() {
        Map<String, Set<String>> map = Maps.newHashMap();

        Set<String> arSet = Sets.newHashSet("V7_splice");
        map.put("AR", arSet);

        Set<String> ntrk1Set = Sets.newHashSet("V3_splice");
        map.put("NTRK1", ntrk1Set);

        return map;
    }
}
