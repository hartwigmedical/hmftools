package com.hartwig.hmftools.vicc.annotation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.common.serve.classification.ImmutableEventClassifierConfig;

import org.jetbrains.annotations.NotNull;

public final class ViccClassificationConfig {

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
    private static final Set<String> HPV_POSITIVE_EVENTS = hpvPositiveEvents();
    private static final Set<String> EBV_POSITIVE_EVENTS = ebvPositiveEvents();
    private static final Map<String, Set<String>> COMBINED_EVENTS_PER_GENE = combinedEventsPerGene();
    private static final Map<String, Set<String>> COMPLEX_EVENTS_PER_GENE = complexEventsPerGene();

    private ViccClassificationConfig() {
    }

    @NotNull
    public static EventClassifierConfig build() {
        return ImmutableEventClassifierConfig.builder()
                .proteinAnnotationExtractor(new ViccProteinAnnotationExtractor())
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
        set.add("EXON");
        set.add("Exon");
        return set;
    }

    @NotNull
    private static Set<String> exonKeywords() {
        Set<String> set = Sets.newHashSet();
        set.add("DELETION");
        set.add("deletion");
        set.add("deletions");
        set.add("deletion/insertion");
        set.add("INSERTION");
        set.add("insertion");
        set.add("insertions");
        set.add("insertions/deletions");
        set.add("mutations");
        set.add("mutation");
        set.add("MUTATION");
        set.add("FRAMESHIFT");
        set.add("proximal");
        return set;
    }

    @NotNull
    private static Set<String> exonBlacklistKeyPhrases() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> specificExonEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("RARE EX 18-21 MUT");
        return set;
    }

    @NotNull
    private static Map<String, Set<String>> fusionPairAndExonsPerGene() {
        Map<String, Set<String>> map = Maps.newHashMap();

        map.put("KIT", Sets.newHashSet("EXON 11 MUTATION", "Exon 11 mutations", "Exon 11 deletions"));
        map.put("MET", Sets.newHashSet("EXON 14 SKIPPING MUTATION"));

        return map;
    }

    @NotNull
    private static Set<String> geneLevelBlacklistKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("exon");
        set.add("EXON");
        set.add("Exon");
        return set;
    }

    @NotNull
    private static Set<String> genericGeneLevelKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("MUTATION");
        set.add("mutant");
        set.add("mut");
        set.add("TRUNCATING MUTATION");
        set.add("Truncating Mutations");
        set.add("feature_truncation");
        set.add("FRAMESHIFT TRUNCATION");
        set.add("FRAMESHIFT MUTATION");
        set.add("ALTERATION");
        set.add("Oncogenic Mutations");
        set.add("oncogenic mutation");
        return set;
    }

    @NotNull
    private static Set<String> activatingGeneLevelKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("Gain-of-function Mutations");
        set.add("Gain-of-Function");
        set.add("act mut");
        set.add("ACTIVATING MUTATION");
        set.add("pos");
        set.add("positive");
        return set;
    }

    @NotNull
    private static Set<String> inactivatingGeneLevelKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("inact mut");
        set.add("biallelic inactivation");
        set.add("Loss Of Function Variant");
        set.add("Loss Of Heterozygosity");
        set.add("DELETERIOUS MUTATION");
        set.add("negative");
        set.add("BIALLELIC INACTIVATION");
        set.add("LOSS-OF-FUNCTION");
        set.add("INACTIVATING MUTATION");
        return set;
    }

    @NotNull
    private static Set<String> amplificationKeywords() {
        Set<String> set = Sets.newHashSet();
        set.add("AMPLIFICATION");
        set.add("Amplification");
        set.add("amplification");
        set.add("amp");
        set.add("overexpression");
        set.add("OVEREXPRESSION");
        set.add("Overexpression");
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
        set.add("EXON");
        set.add("Exon");
        set.add("Ex19");
        set.add("inframe");
        return set;
    }

    @NotNull
    private static Set<String> deletionKeywords() {
        Set<String> set = Sets.newHashSet();
        set.add("Deletion");
        set.add("deletion");
        set.add("DELETION");
        set.add("del");
        set.add("undexpression");
        set.add("UNDEREXPRESSION");
        set.add("loss");
        set.add("LOSS");
        return set;
    }

    @NotNull
    private static Set<String> deletionKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("dec exp");
        set.add("Copy Number Loss");
        return set;
    }

    @NotNull
    private static Set<String> exonicDelDupFusionKeyPhrases() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> exonicDelDupFusionEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("EGFRvII");
        set.add("EGFRvIII");
        set.add("EGFRvV");
        set.add("EGFR-KDD");
        return set;
    }

    @NotNull
    private static Set<String> fusionPairEventsToSkip() {
        Set<String> set = Sets.newHashSet();
        set.add("AR-V7");
        set.add("Gain-of-Function");
        set.add("LOSS-OF-FUNCTION");
        set.add("LCS6-variant");
        set.add("DI842-843VM");
        set.add("FLT3-ITD");
        return set;
    }

    @NotNull
    private static Set<String> promiscuousFusionKeyPhrases() {
        Set<String> set = Sets.newHashSet();
        set.add("Fusion");
        set.add("fusion");
        set.add("FUSION");
        set.add("Fusions");
        set.add("FUSIONS");
        set.add("REARRANGEMENT");
        set.add("rearrange");
        return set;
    }

    @NotNull
    private static Set<String> microsatelliteUnstableEvents() {
        Set<String> set = Sets.newHashSet();
        set.add("Microsatellite Instability-High");
        return set;
    }

    @NotNull
    private static Set<String> microsatelliteStableEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> highTumorMutationalLoadEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> lowTumorMutationalLoadEvents() {
        return Sets.newHashSet();
    }

    @NotNull
    private static Set<String> hrDeficiencyEvents() {
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
        Map<String, Set<String>> map = Maps.newHashMap();

        map.put("EGFR", Sets.newHashSet("Ex19 del L858R"));
        map.put("BRAF", Sets.newHashSet("p61BRAF-V600E", "V600E AMPLIFICATION"));

        return map;
    }

    @NotNull
    private static Map<String, Set<String>> complexEventsPerGene() {
        Map<String, Set<String>> map = Maps.newHashMap();

        Set<String> alkSet = Sets.newHashSet("ALK inframe insertion (1151T)");
        map.put("ALK", alkSet);

        Set<String> arSet = Sets.newHashSet("ARv567es", "SPLICE VARIANT 7", "AR-V7");
        map.put("AR", arSet);

        Set<String> brafSet = Sets.newHashSet("DEL 485-490", "L485_P490>Y", "BRAF p.L485_P490");
        map.put("BRAF", brafSet);

        Set<String> brca1Set = Sets.newHashSet("Alu insertion");
        map.put("BRCA1", brca1Set);

        Set<String> casp8Set = Sets.newHashSet("CASP8L");
        map.put("CASP8", casp8Set);

        Set<String> cebpaSet = Sets.newHashSet("N-TERMINAL FRAME SHIFT");
        map.put("CEBPA", cebpaSet);

        Set<String> ccnd1Set = Sets.newHashSet("256_286trunc");
        map.put("CCND1", ccnd1Set);

        Set<String> ccnd3Set = Sets.newHashSet("D286_L292trunc");
        map.put("CCND3", ccnd3Set);

        Set<String> chek2Set = Sets.newHashSet("1100DELC", "IVS2+1G>A");
        map.put("CHEK2", chek2Set);

        Set<String> dnmt3bSet = Sets.newHashSet("DNMT3B7");
        map.put("DNMT3B", dnmt3bSet);

        Set<String> dpydSet = Sets.newHashSet("DPYD splice donor variant", "DPYD*13 HOMOZYGOSITY", "DPYD*2A HOMOZYGOSITY");
        map.put("DPYD", dpydSet);

        Set<String> egfrSet = Sets.newHashSet("EGFR inframe deletion (30-336)",
                "EGFR inframe insertion (769-770)",
                "EGFR inframe deletion (6-273)",
                "EGFR CTD");
        map.put("EGFR", egfrSet);

        Set<String> eif1axSet = Sets.newHashSet("A113_splice");
        map.put("EIF1AX", eif1axSet);

        Set<String> epcamSet = Sets.newHashSet("3' EXON DELETION");
        map.put("EPCAM", epcamSet);

        Set<String> erbb2Set = Sets.newHashSet("P780INS", "DEL 755-759", "KINASE DOMAIN MUTATION");
        map.put("ERBB2", erbb2Set);

        Set<String> ezh2Set = Sets.newHashSet("INTRON 6 MUTATION");
        map.put("EZH2", ezh2Set);

        Set<String> fli1Set = Sets.newHashSet("EWSR1-FLI1 Type 1");
        map.put("FLI1", fli1Set);

        Set<String> flt3Set = Sets.newHashSet("FLT3-ITD", "ITD", "FLT3 internal tandem duplications", "TKD MUTATION");
        map.put("FLT3", flt3Set);

        Set<String> hlaaSet = Sets.newHashSet("596_619splice");
        map.put("HLA-A", hlaaSet);

        Set<String> jak2Set = Sets.newHashSet("F547 SPLICE SITE MUTATION");
        map.put("JAK2", jak2Set);

        Set<String> kitSet =
                Sets.newHashSet("INTERNAL DUPLICATION", "3' UTR MUTATION", "KIT inframe deletion (577-579)", "KIT inframe deletion (V560)");
        map.put("KIT", kitSet);

        Set<String> krasSet = Sets.newHashSet("LCS6-variant");
        map.put("KRAS", krasSet);

        Set<String> map2k1Set = Sets.newHashSet("MAP2K1 inframe deletion (56-60)");
        map.put("MAP2K1", map2k1Set);

        Set<String> metSet = Sets.newHashSet("MET kinase domain mutation",
                "963_D1010splice",
                "X963_splice",
                "981_1028splice",
                "X1006_splice",
                "X1007_splice",
                "X1008_splice",
                "X1009_splice");
        map.put("MET", metSet);

        Set<String> mlh1Set = Sets.newHashSet("C.790+1G>A");
        map.put("MLH1", mlh1Set);

        Set<String> notch1Set = Sets.newHashSet("NOTCH1 activating mutation in Cterm-PEST domain",
                "NOTCH2 activating mutation (missense in TAD or truncating in Cterm-PEST domain)",
                "Truncating Mutations in the PEST Domain",
                "Truncating Mutations Upstream of Transactivation Domain");
        map.put("NOTCH1", notch1Set);

        Set<String> notch2Set = Sets.newHashSet("2010_2471trunc",
                "1_2009trunc",
                "NOTCH2 activating mutation (missense in TAD or truncating in Cterm-PEST domain)");
        map.put("NOTCH2", notch2Set);

        Set<String> ntrk1Set = Sets.newHashSet("TRKAIII Splice Variant");
        map.put("NTRK1", ntrk1Set);

        Set<String> pdgfraSet = Sets.newHashSet("PDGFRA inframe deletion (I843)");
        map.put("PDGFRA", pdgfraSet);

        Set<String> pik3r1Set = Sets.newHashSet("X582_splice", "X475_splice", "X434_splice");
        map.put("PIK3R1", pik3r1Set);

        Set<String> pmlSet = Sets.newHashSet("B2 DOMAIN MUTATION");
        map.put("PML", pmlSet);

        Set<String> poleSet = Sets.newHashSet("POLE (268-471)");
        map.put("POLE", poleSet);

        Set<String> ppm1dSet = Sets.newHashSet("422_605trunc");
        map.put("PPM1D", ppm1dSet);

        Set<String> ptch1Set = Sets.newHashSet("LOH");
        map.put("PTCH1", ptch1Set);

        Set<String> runx1Set = Sets.newHashSet("R135FSX177", "T148HFSX9");
        map.put("RUNX1", runx1Set);

        Set<String> tertSet = Sets.newHashSet("TERT promoters core", "Promoter Mutations", "PROMOTER MUTATION");
        map.put("TERT", tertSet);

        Set<String> tgfbr1Set = Sets.newHashSet("TGFBR1*6A");
        map.put("TGFBR1", tgfbr1Set);

        Set<String> tp53Set = Sets.newHashSet("DNA binding domain deletions",
                "DNA binding domain insertions",
                "DNA binding domain missense mutations",
                "DNA BINDING DOMAIN MUTATION");
        map.put("TP53", tp53Set);

        Set<String> tpmtSet = Sets.newHashSet("TPMT splice acceptor variant");
        map.put("TPMT", tpmtSet);

        Set<String> tymsSet = Sets.newHashSet("5' TANDEM REPEAT");
        map.put("TYMS", tymsSet);

        Set<String> ugt1a1Set = Sets.newHashSet("UGT1A1*28", "UGT1A1*60");
        map.put("UGT1A1", ugt1a1Set);

        Set<String> vhlSet = Sets.newHashSet("R108ins (c.324InsCGC)",
                "3'UTR alteration (c.639+10C>G)",
                "3'UTR alteration (c.642+70C>A)",
                "Splicing alteration (c.341-2A>C)",
                "Splicing alteration (c.463+1G>C)",
                "Splicing alteration (c.463+2C>T)",
                "Splicing alteration (c.464-2A>G)",
                "Splicing alteration (c.464-2A>T)",
                "Splicing alteration (c.464-1G>A)",
                "Splicing alteration (c.464-1G>C)");
        map.put("VHL", vhlSet);

        return map;
    }
}
