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
    private static final Set<String> EXONIC_DEL_DUP_FUSION_EVENTS = exonicDelDupFusionEvents();
    private static final Set<String> FUSION_PAIR_EVENTS_TO_SKIP = fusionPairEventsToSkip();
    private static final Set<String> PROMISCUOUS_FUSION_KEY_PHRASES = promiscuousFusionKeyPhrases();
    private static final Set<String> MICROSATELLITE_UNSTABLE_EVENTS = microsatelliteUnstableEvents();
    private static final Set<String> MICROSATELLITE_STABLE_EVENTS = microsatelliteStableEvents();
    private static final Set<String> HIGH_TUMOR_MUTATIONAL_LOAD_EVENTS = highTumorMutationalLoadEvents();
    private static final Set<String> LOW_TUMOR_MUTATIONAL_LOAD_EVENTS = lowTumorMutationalLoadEvents();
    private static final Set<String> HR_DEFICIENCY_EVENTS = hrDeficiencyEvents();
    private static final Set<String> HPV_POSITIVE_EVENTS = hpvPositiveEvents();
    private static final Set<String> EBV_POSITIVE_EVENTS = ebvPositiveEvents();
    private static final Map<String, Set<String>> COMBINED_EVENTS_PER_GENE = combinedEventsPerGene();
    private static final Map<String, Set<String>> COMPLEX_EVENTS_PER_GENE = complexEventsPerGene();

    private CkbClassificationConfig() {
    }

    @NotNull
    public static EventClassifierConfig build() {
        return ImmutableEventClassifierConfig.builder()
                .proteinAnnotationExtractor(new ProteinAnnotationExtractor())
                .exonIdentifiers(EXON_IDENTIFIERS)
                .exonKeywords(EXON_KEYWORDS)
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
                .exonicDelDupFusionEvents(EXONIC_DEL_DUP_FUSION_EVENTS)
                .fusionPairEventsToSkip(FUSION_PAIR_EVENTS_TO_SKIP)
                .promiscuousFusionKeyPhrases(PROMISCUOUS_FUSION_KEY_PHRASES)
                .microsatelliteUnstableEvents(MICROSATELLITE_UNSTABLE_EVENTS)
                .microsatelliteStableEvents(MICROSATELLITE_STABLE_EVENTS)
                .highTumorMutationalLoadEvents(HIGH_TUMOR_MUTATIONAL_LOAD_EVENTS)
                .lowTumorMutationalLoadEvents(LOW_TUMOR_MUTATIONAL_LOAD_EVENTS)
                .hrDeficiencyEvents(HR_DEFICIENCY_EVENTS)
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
    private static Set<String> specificExonEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Map<String, Set<String>> fusionPairAndExonsPerGene() {
        Map<String, Set<String>> map = Maps.newHashMap();
        map.put("KIT", Sets.newHashSet("exon  11 del", "exon 11"));
        map.put("MET", Sets.newHashSet("del exon 14", "exon 14"));

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
        set.add("positive");
        return set;
    }

    @NotNull
    private static Set<String> inactivatingGeneLevelKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("inact mut");
        set.add("negative");
        set.add("LOH");
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
    private static Set<String> exonicDelDupFusionEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> fusionPairEventsToSkip() {
        Set<String> set = Sets.newHashSet();
        return set;
    }

    @NotNull
    private static Set<String> promiscuousFusionKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("fusion promisuous");
        set.add("rearrange");
        return set;
    }

    @NotNull
    private static Set<String> microsatelliteUnstableEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("MSI HIGH");
        return set;
    }

    @NotNull
    private static Set<String> microsatelliteStableEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("MSI LOW");
        return set;
    }

    @NotNull
    private static Set<String> highTumorMutationalLoadEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("TumMutLoad HIGH");
        return set;
    }

    @NotNull
    private static Set<String> lowTumorMutationalLoadEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("TumMutLoad LOW");
        return set;
    }

    @NotNull
    private static Set<String> hrDeficiencyEvents() {
        Set<String> set = Sets.newHashSet();
        return set;
    }

    @NotNull
    private static Set<String> hpvPositiveEvents() {
        Set<String> set = Sets.newHashSet();
        return set;
    }

    @NotNull
    private static Set<String> ebvPositiveEvents() {
        Set<String> set = Sets.newHashSet();
        return set;
    }

    @NotNull
    private static Map<String, Set<String>> combinedEventsPerGene() {
        Map<String, Set<String>> map = Maps.newHashMap();
        return map;
    }

    @NotNull
    private static Map<String, Set<String>> complexEventsPerGene() {
        Map<String, Set<String>> map = Maps.newHashMap();

        Set<String> arSet = Sets.newHashSet("V7_splice");
        map.put("AR", arSet);

        Set<String> atmSet = Sets.newHashSet("G2718_K2756del", "N2326_K2363del", "R2506_N2543del");
        map.put("ATM", atmSet);

        Set<String> axin2Set = Sets.newHashSet("E745_S762del", "S217_E823del");
        map.put("AXIN2", axin2Set);

        Set<String> brafSet = Sets.newHashSet("A81_D380del",
                "A81_M438del",
                "G203_G393del",
                "V169_D380del",
                "V169_G327del",
                "V47_D380del",
                "V47_G327del",
                "V47_G393del",
                "V47_M438del");
        map.put("BRAF", brafSet);

        Set<String> brca1Set = Sets.newHashSet("E427_S713del", "S377_N417del");
        map.put("BRCA1", brca1Set);

        Set<String> brca2Set = Sets.newHashSet("S917_H1918del");
        map.put("BRCA2", brca2Set);

        Set<String> brd4Set = Sets.newHashSet("Q762_P780del");
        map.put("BRD4", brd4Set);

        Set<String> btkSet = Sets.newHashSet("P204_Y263del");
        map.put("BTK", btkSet);

        Set<String> cblSet = Sets.newHashSet("E366_K382del");
        map.put("CBL", cblSet);

        Set<String> chek2Set = Sets.newHashSet("D265_H282del", "E107_K197del");
        map.put("CHEK2", chek2Set);

        Set<String> cic2Set = Sets.newHashSet("R1464_M1519del", "W24_V311del");
        map.put("CIC", cic2Set);

        Set<String> ctnna1Set = Sets.newHashSet("V36_Q196del");
        map.put("CTNNA1", ctnna1Set);

        Set<String> ctnnb1Set = Sets.newHashSet("A21_A149del",
                "A21_A152del",
                "A5_A80del",
                "A5_Q143del",
                "H24_Y142del",
                "L10_N141del",
                "M12_D144del",
                "M8_L132del",
                "P16_K133del",
                "T3_A126del",
                "V22_A97del",
                "V22_D145del",
                "V22_G38del",
                "V22_Y64del");
        map.put("CTNNB1", ctnnb1Set);

        Set<String> egfrSet = Sets.newHashSet("G696_P1033dup", "G983_G1054del", "T34_A289del", "V1010_D1152del", "V30_R297delinsG");
        map.put("EGFR", egfrSet);

        Set<String> erbb3Set = Sets.newHashSet("E928fs*16");
        map.put("ERBB3", erbb3Set);

        Set<String> fancfSet = Sets.newHashSet("L285_V302del");
        map.put("FANCF", fancfSet);

        Set<String> frs2Set = Sets.newHashSet("F251_Q270del", "P232_G272del", "P232_R295del");
        map.put("FRS2", frs2Set);

        Set<String> gata1Set = Sets.newHashSet("D65_C228del");
        map.put("GATA1", gata1Set);

        Set<String> jak1Set = Sets.newHashSet("E966_K989del");
        map.put("JAK1", jak1Set);

        Set<String> kitSet = Sets.newHashSet("E554_D572delinsA",
                "M552_D572del",
                "V555_I571del",
                "V555_L576del",
                "V555_P573del",
                "V556_L576del",
                "V556_T574del",
                "V560_L576del",
                "V560_Y578del");
        map.put("KIT", kitSet);

        Set<String> lrp1bSet = Sets.newHashSet("D155_I197del");
        map.put("LRP1B", lrp1bSet);

        Set<String> mlh1Set = Sets.newHashSet("E633_E663del", "P578_E632del");
        map.put("MLH1", mlh1Set);

        Set<String> notch1Set = Sets.newHashSet("S2499_F2554del", "V1110_S1723del");
        map.put("NOTCH1", notch1Set);

        Set<String> ntrk1Set = Sets.newHashSet("V3_splice");
        map.put("NTRK1", ntrk1Set);

        Set<String> pdgfraSet = Sets.newHashSet("C456_R481del", "Y375_K455del");
        map.put("PDGFRA", pdgfraSet);

        Set<String> pik3r1Set = Sets.newHashSet("I442_Y462del");
        map.put("PIK3R1", pik3r1Set);

        Set<String> plcg2Set = Sets.newHashSet("A686_W806del", "W646_R685del");
        map.put("PLCG2", plcg2Set);

        Set<String> ptk2Set = Sets.newHashSet("E956_E982del");
        map.put("PTK2", ptk2Set);

        Set<String> rb1Set = Sets.newHashSet("T738_R775del");
        map.put("RB1", rb1Set);

        Set<String> stk11Set = Sets.newHashSet("E98_G155del");
        map.put("STK11", stk11Set);

        Set<String> znf703Set = Sets.newHashSet("K295_S327del");
        map.put("ZNF703", znf703Set);

        return map;
    }
}
