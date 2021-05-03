package com.hartwig.hmftools.serve.sources.vicc.curation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class FeatureCurationFactory {

    static final Map<FeatureCurationKey, FeatureCurationValues> FEATURE_MAPPINGS = Maps.newHashMap();

    static final Set<FeatureCurationKey> FEATURE_BLACKLIST = Sets.newHashSet();

    private FeatureCurationFactory() {
    }

    static {
        populateCGICuration();
        populateCIViCCuration();
        populateJaxCuration();
        populateOncoKBCuration();
    }

    private static void populateCGICuration() {
        // These are not wrong but lead to inconsistencies downstream
        FEATURE_MAPPINGS.put(cgi("PIK3R1", null, "PIK3R1 p.E439delE"), curation("PIK3R1", "PIK3R1 p.E439del"));
        FEATURE_MAPPINGS.put(cgi("PIK3R1", null, "PIK3R1 p.D560_S565delDKRMNS"), curation("PIK3R1", "PIK3R1 p.D560_S565del"));
        FEATURE_MAPPINGS.put(cgi("PIK3R1", null, "PIK3R1 p.T576delT"), curation("PIK3R1", "PIK3R1 p.T576del"));

        // These fusions need to be swapped around.
        FEATURE_MAPPINGS.put(cgi("ABL1", null, "ABL1-BCR fusion"), curation("BCR", "BCR-ABL1 fusion"));
        FEATURE_MAPPINGS.put(cgi("BCR", null, "ABL1-BCR fusion"), curation("BCR", "BCR-ABL1 fusion"));
        FEATURE_MAPPINGS.put(cgi("COL1A1", null, "PDGFB-COL1A1 fusion"), curation("COL1A1", "COL1A1-PDGFB fusion"));
        FEATURE_MAPPINGS.put(cgi("PDGFB", null, "PDGFB-COL1A1 fusion"), curation("COL1A1", "COL1A1-PDGFB fusion"));
        FEATURE_MAPPINGS.put(cgi("FIP1L1", null, "PDGFRA-FIP1L1 fusion"), curation("FIP1L1", "FIP1L1-PDGFRA fusion"));
        FEATURE_MAPPINGS.put(cgi("PDGFRA", null, "PDGFRA-FIP1L1 fusion"), curation("FIP1L1", "FIP1L1-PDGFRA fusion"));
        FEATURE_MAPPINGS.put(cgi("PRKCH", null, "PRKCH amplification + ABL1-BCR fusion"),
                curation("PRKCH", "PRKCH amplification + BCR-ABL1 fusion"));

        // These fusions are on synonym genes
        FEATURE_MAPPINGS.put(cgi("BRD4", null, "BRD4-C15orf55 fusion"), curation("BRD4", "BRD4-NUTM1 fusion"));
        FEATURE_MAPPINGS.put(cgi("C15orf55", null, "BRD4-C15orf55 fusion"), curation("BRD4", "BRD4-NUTM1 fusion"));
        FEATURE_MAPPINGS.put(cgi("MLL", null, "MLL fusion"), curation("KMT2A", "KMT2A fusion"));

        // These genes are on synonym genes
        FEATURE_MAPPINGS.put(cgi("MLL2", null, "MLL2 oncogenic mutation"), curation("KMT2D", "KMT2D oncogenic mutation"));

        // Improve for consistency
        FEATURE_MAPPINGS.put(cgi("MET", null, "MET (Y1230C;Y1235D)"), curation("MET", "MET (Y1230C,Y1235D)"));

        // Variants that probably exist on another transcript
        FEATURE_MAPPINGS.put(cgi("FGFR3", null, "FGFR3 (K650)"), curation("FGFR3", "FGFR3 (K652)"));
        FEATURE_MAPPINGS.put(cgi("GNAS", null, "GNAS (R201)"), curation("GNAS", "GNAS (R844)"));
    }

    private static void populateCIViCCuration() {
        // Protein annotation is not formally correct, update to make it correct.
        FEATURE_MAPPINGS.put(civic("EGFR", "ENST00000275493", "V769_770insASV"), curation("EGFR", "V769_D770insASV"));
        FEATURE_MAPPINGS.put(civic("ERBB2", "ENST00000269571", "M774INSAYVM"), curation("ERBB2", "M774_A775INSAYVM"));

        // Fusion needs to be flipped around
        FEATURE_MAPPINGS.put(civic("BRAF", "ENST00000288602", "BRAF-CUL1"), curation("BRAF", "CUL1-BRAF"));

        // Fusions where gene is a synonym of our gene
        FEATURE_MAPPINGS.put(civic("FGFR1", null, "ZNF198-FGFR1"), curation("ZMYM2", "ZMYM2-FGFR1"));
        FEATURE_MAPPINGS.put(civic("KMT2A", null, "MLL-MLLT3"), curation("KMT2A", "KMT2A-MLLT3"));

        // Map specific EGFR fusion do a different name to be consistent with other sources.
        FEATURE_MAPPINGS.put(civic("EGFR", "ENST00000275493", "VIII"), curation("EGFR", "EGFRvIII"));

        // These genes are on synonym genes
        FEATURE_MAPPINGS.put(civic("MRE11", "ENST00000323929", "LOSS"), curation("MRE11A", "LOSS"));
        FEATURE_MAPPINGS.put(civic("MRE11", "ENST00000323929", "FRAMESHIFT MUTATION"), curation("MRE11A", "FRAMESHIFT MUTATION"));

        // Fusions are not correct
        FEATURE_MAPPINGS.put(civic("ABL1", "ENST00000318560", "BCR-ABL G398R"), curation("BCR", "BCR-ABL1 G398R"));
        FEATURE_MAPPINGS.put(civic("ABL1", "ENST00000372348", "BCR-ABL T334I"), curation("BCR", "BCR-ABL1 T334I"));
        FEATURE_MAPPINGS.put(civic("ABL1", "ENST00000318560", "BCR-ABL E255K"), curation("BCR", "BCR-ABL1 E255K"));
        FEATURE_MAPPINGS.put(civic("ABL1", "ENST00000318560", "BCR-ABL F317L"), curation("BCR", "BCR-ABL1 F317L"));
        FEATURE_MAPPINGS.put(civic("ABL1", "ENST00000318560", "BCR-ABL F486S"), curation("BCR", "BCR-ABL1 F486S"));
        FEATURE_MAPPINGS.put(civic("ABL1", "ENST00000305877", "BCR-ABL"), curation("BCR", "BCR-ABL1"));

        // Variants that don't exist
        FEATURE_BLACKLIST.add(civic("NOTCH1", "ENST00000277541", "S2275FS"));
        FEATURE_BLACKLIST.add(civic("NOTCH1", "ENST00000277541", "V2444FS"));
        FEATURE_BLACKLIST.add(civic("VHL", null, "L132fs (c.395delA)"));
        FEATURE_BLACKLIST.add(civic("GNAS", "ENST00000371100", "T393C"));
        FEATURE_BLACKLIST.add(civic("TERT", "ENST00000310581", "C228T"));
        FEATURE_BLACKLIST.add(civic("VHL", null, "P81delRVV (c.243_251delGCGCGTCGT)"));
        FEATURE_BLACKLIST.add(civic("VHL", null, "A121I (c.364_365GC>AT)"));
        FEATURE_BLACKLIST.add(civic("VHL", null, "N167T (c.392A>C)"));
        FEATURE_BLACKLIST.add(civic("VHL", null, "X214L (c.641G>T)"));

        // Variant only possible through an MNV which spans an intron, and inconsistent with its coding impact.
        FEATURE_BLACKLIST.add(civic("VHL", null, "V155C (c.462delA)"));

        // Variant unlikely to be real as it spans multiple exons
        FEATURE_BLACKLIST.add(civic("VHL", null, "G114dup (c.342dupGGT)"));

        // Variant hard to interpret as it crossing exonic boundary (in protein space)
        FEATURE_BLACKLIST.add(civic("VHL", null, "G114fs (c.341delCGTTTCCAACAATTTCTCGGTGT)"));

        // Variant hard to interpret in protein space as it changes an amino acid followed by a frameshift.
        FEATURE_BLACKLIST.add(civic("VHL", "ENST00000256474", "M54IFS (c.162_166delGGAGG)"));

        // Variant has an inconsistent formatting. Not sure how to correct
        FEATURE_BLACKLIST.add(civic("PDGFRA", "ENST00000257290", "DI842-843VM"));
    }

    private static void populateJaxCuration() {
        // Update protein annotation to be correct (should be capitalized).
        FEATURE_MAPPINGS.put(jax("PIK3CA", null, "PIK3CA E545k "), curation("PIK3CA", "PIK3CA E545K "));

        // These mappings are to work around the missing transcripts in JAX.
        // We map every mutation that appears on both the canonical + non-canonical form to its canonical form in ensembl.
        FEATURE_MAPPINGS.put(jax("EZH2", null, "EZH2 Y641F "), curation("EZH2", "EZH2 Y646F "));
        FEATURE_MAPPINGS.put(jax("EZH2", null, "EZH2 Y641H "), curation("EZH2", "EZH2 Y646H "));
        FEATURE_MAPPINGS.put(jax("EZH2", null, "EZH2 Y641N "), curation("EZH2", "EZH2 Y646N "));
        FEATURE_MAPPINGS.put(jax("EZH2", null, "EZH2 Y641S "), curation("EZH2", "EZH2 Y646S "));
        FEATURE_MAPPINGS.put(jax("FGFR2", null, "FGFR2 V564I "), curation("FGFR2", "FGFR2 V565I "));
        FEATURE_MAPPINGS.put(jax("FGFR3", null, "FGFR3 Y373C "), curation("FGFR3", "FGFR3 Y375C "));
        FEATURE_MAPPINGS.put(jax("FGFR3", null, "FGFR3 K650E "), curation("FGFR3", "FGFR3 K652E "));
        FEATURE_MAPPINGS.put(jax("MET", null, "MET L1195V "), curation("MET", "MET L1213V "));
        FEATURE_MAPPINGS.put(jax("MET", null, "MET Y1230H "), curation("MET", "MET Y1248H "));
        FEATURE_MAPPINGS.put(jax("MET", null, "MET M1250T "), curation("MET", "MET M1268T "));

        // These mappings are identical but used concurrently. Confirmed to be replaced by I740_K745dup
        FEATURE_MAPPINGS.put(jax("EGFR", null, "EGFR I744_K745insKIPVAI "), curation("EGFR", "EGFR K745_E746insIPVAIK "));
        FEATURE_MAPPINGS.put(jax("KIT", null, "KIT V559del "), curation("KIT", "KIT V560del "));

        // The below variants in FLT3 are from a paper where an additional R was added in the ref sequence, shifting all AAs by one position.
        // This has been corrected in current live CKB.
        FEATURE_BLACKLIST.add(jax("FLT3", null, "FLT3 L611_E612insCSSDNEYFYVDFREYEYDLKWEFPRENL "));
        FEATURE_BLACKLIST.add(jax("FLT3", null, "FLT3 E612_F613insGYVDFREYEYDLKWEFRPRENLEF "));

        // The transcript that should have this mutation (ENST00000507379) is annotated as 3' truncated with only 1135 AAs in ensembl)
        FEATURE_BLACKLIST.add(jax("APC", null, "APC S1197* "));

        // The below is pending investigation by JAX, possibly a mistake by the paper.
        FEATURE_BLACKLIST.add(jax("PIK3CA", null, "PIK3CA R425L "));

        // Variant that doesn't seem to exist.
        FEATURE_BLACKLIST.add(jax("PTEN", null, "PTEN Y86fs "));

        // Variants spanning multiple exons in our interpretation
        FEATURE_BLACKLIST.add(jax("PTEN", null, "PTEN R55fs*1 "));
        FEATURE_BLACKLIST.add(jax("PTEN", null, "PTEN Y27fs*1 "));

        // Variant hard to interpret as it crossing exonic boundary
        FEATURE_BLACKLIST.add(jax("VHL", null, "VHL V155fs "));
    }

    private static void populateOncoKBCuration() {
        // The dashes in these keys are not the usual dashes!
        FEATURE_MAPPINGS.put(oncoKb("FGFR1", "ENST00000425967", "ERLIN2–FGFR1 Fusion"), curation("ERLIN2", "ERLIN2-FGFR1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("FGFR2", "ENST00000358487", "FGFR2–PPHLN1 Fusion"), curation("FGFR2", "FGFR2-PPHLN1 Fusion"));

        // The spaces are not consistent
        FEATURE_MAPPINGS.put(oncoKb("FGFR3", "ENST00000260795", "FGFR3 - BAIAP2L1 Fusion"), curation("FGFR3", "FGFR3-BAIAP2L1 Fusion"));

        // These fusions need to be flipped.
        FEATURE_MAPPINGS.put(oncoKb("CCND1", "ENST00000227507", "CCND1-IGH Fusion"), curation("IGH", "IGH-CCND1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("ROS1", "ENST00000368508", "ROS1-CD74 Fusion"), curation("CD74", "CD74-ROS1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("NSD1", "ENST00000439151", "NSD1-NUP98 Fusion"), curation("NUP98", "NUP98-NSD1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("RET", "ENST00000355710", "RET-CCDC6 Fusion"), curation("CCDC6", "CCDC6-RET Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("EP300", "ENST00000263253", "EP300-MOZ Fusion"), curation("KAT6A", "KAT6A-EP300 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("EP300", "ENST00000263253", "EP300-MLL Fusion"), curation("KMT2A", "KMT2A-EP300 Fusion"));

        // These fusions are incorrect but recoverable
        FEATURE_MAPPINGS.put(oncoKb("PAX8", "ENST00000263334", "PAX8-PPARγ Fusion"), curation("PAX8", "PAX8-PPARG Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("ALK", "ENST00000389048", "NPM-ALK Fusion"), curation("NPM1", "NPM1-ALK Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("BRD4", "ENST00000263377", "BRD4-NUT Fusion"), curation("BRD4", "BRD4-NUTM1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("NOTCH1", "ENST00000277541", "SEC16A1-NOTCH1 Fusion"), curation("SEC16A", "SEC16A-NOTCH1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("FGFR1", "ENST00000425967", "FGFR1OP1-FGFR1 Fusion"), curation("FGFR1OP", "FGFR1OP-FGFR1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("NKX2-1", "ENST00000354822", "IGH-NKX2 Fusion"), curation("NKX2-1", "IGH-NKX2-1 Fusion"));

        // These fusions are on synonym genes
        FEATURE_MAPPINGS.put(oncoKb("FGFR1", "ENST00000425967", "ZNF198-FGFR1 Fusion"), curation("ZMYM2", "ZMYM2-FGFR1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("TET1", "ENST00000373644", "MLL-TET1 Fusion"), curation("KMT2A", "KMT2A-TET1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("ROS1", "ENST00000368508", "FIG-ROS1 Fusion"), curation("GOPC", "GOPC-ROS1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("PDGFRB", "ENST00000261799", "GPIAP1-PDGFRB Fusion"), curation("CAPRIN1", "CAPRIN1-PDGFRB Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("PDGFRB", "ENST00000261799", "KIAA1509-PDGFRB Fusion"), curation("CCDC88C", "CCDC88C-PDGFRB Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("JAK2", "ENST00000381652", "TEL-JAK2 Fusion"), curation("ETV6", "ETV6-JAK2 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("RUNX1", "ENST00000300305", "TEL-RUNX1 Fusion"), curation("ETV6", "ETV6-RUNX1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("RUNX1", "ENST00000300305", "RUNX1-EVI1 Fusion"), curation("RUNX1", "RUNX1-MECOM Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("FGFR1", "ENST00000425967", "CEP110-FGFR1 Fusion"), curation("CNTRL", "CNTRL-FGFR1 Fusion"));
        FEATURE_MAPPINGS.put(oncoKb("FGFR2", "ENST00000358487", "FGFR2-KIAA1967 Fusion"), curation("FGFR2", "FGFR2-CCAR2 Fusion"));

        // These are inconsistent or improperly aligned variants.
        FEATURE_MAPPINGS.put(oncoKb("BRAF", "ENST00000288602", "T599insTT"), curation("BRAF", "T599_V600insTT"));
        FEATURE_MAPPINGS.put(oncoKb("EGFR", "ENST00000275493", "E746_T751insIP"), curation("EGFR", "E746_L747insIP"));
        FEATURE_MAPPINGS.put(oncoKb("EGFR", "ENST00000275493", "H773insLGNP"), curation("EGFR", "H773_V774insLGNP"));
        FEATURE_MAPPINGS.put(oncoKb("EPAS1", "ENST00000263734", "533_534del"), curation("EPAS1", "I533_P534del"));
        FEATURE_MAPPINGS.put(oncoKb("EPAS1", "ENST00000263734", "534_536del"), curation("EPAS1", "P534_D536del"));
        FEATURE_MAPPINGS.put(oncoKb("KIT", "ENST00000288135", "V559del"), curation("KIT", "V560del"));
        FEATURE_MAPPINGS.put(oncoKb("PTEN", "ENST00000371953", "I32del"), curation("PTEN", "I33del"));
        FEATURE_MAPPINGS.put(oncoKb("RIT1", "ENST00000368323", "T76insTLDT"), curation("RIT1", "T76_A77insTLDT"));

        // Not clear what TRA and TRB are
        FEATURE_BLACKLIST.add(oncoKb("NKX2-1", "ENST00000354822", "TRA-NKX2-1 Fusion"));
        FEATURE_BLACKLIST.add(oncoKb("NKX2-1", "ENST00000354822", "TRB-NKX2-1 Fusion"));

        // Not sure what "Delta" means in this context.
        FEATURE_BLACKLIST.add(oncoKb("NTRK1", "ENST00000524377", "Delta-NTRK1 Fusion"));

        // Variants are unlikely as they span multiple exons (and hence are more fusions than inframes)
        FEATURE_BLACKLIST.add(oncoKb("PDGFRA", "ENST00000257290", "E311_K312del"));
        FEATURE_BLACKLIST.add(oncoKb("ETV6", "ENST00000396373", "385_418del"));

        // Variants that we couldn't quite figure out
        FEATURE_BLACKLIST.add(oncoKb("KIT", "ENST00000288135", "T574insTQLPYD"));

        // Variants that don't exist
        //  - Below is called a "silent promoter" according to https://pubmed.ncbi.nlm.nih.gov/11606402/
        FEATURE_BLACKLIST.add(oncoKb("APC", "ENST00000257430", "A290T"));
        FEATURE_BLACKLIST.add(oncoKb("AR", "ENST00000374690", "F876L"));
        FEATURE_BLACKLIST.add(oncoKb("BRCA1", "ENST00000357654", "T47D"));
        FEATURE_BLACKLIST.add(oncoKb("BRCA1", "ENST00000357654", "E60L"));
        FEATURE_BLACKLIST.add(oncoKb("BRCA2", "ENST00000380152", "S1670A"));
        FEATURE_BLACKLIST.add(oncoKb("BRCA2", "ENST00000380152", "R2304C"));
        FEATURE_BLACKLIST.add(oncoKb("BRCA2", "ENST00000380152", "N2829R"));
        FEATURE_BLACKLIST.add(oncoKb("CARD11", "ENST00000396946", "G116S"));
        FEATURE_BLACKLIST.add(oncoKb("CARD11", "ENST00000396946", "F123I"));
        FEATURE_BLACKLIST.add(oncoKb("CARD11", "ENST00000396946", "E127G"));
        FEATURE_BLACKLIST.add(oncoKb("CARM1", "ENST00000327064", "S217A"));
        FEATURE_BLACKLIST.add(oncoKb("CARM1", "ENST00000327064", "S217C"));
        FEATURE_BLACKLIST.add(oncoKb("CARM1", "ENST00000327064", "S217E"));
        FEATURE_BLACKLIST.add(oncoKb("CARM1", "ENST00000327064", "S217T"));
        FEATURE_BLACKLIST.add(oncoKb("CASP8", "ENST00000358485", "C248T"));
        FEATURE_BLACKLIST.add(oncoKb("CCND3", "ENST00000372991", "T286A"));
        FEATURE_BLACKLIST.add(oncoKb("CDK12", "ENST00000447079", "K765R"));
        FEATURE_BLACKLIST.add(oncoKb("CDK12", "ENST00000447079", "D887N"));
        FEATURE_BLACKLIST.add(oncoKb("CDKN2A", "ENST00000304494", "P73L"));
        FEATURE_BLACKLIST.add(oncoKb("CDKN2A", "ENST00000304494", "R79P"));
        FEATURE_BLACKLIST.add(oncoKb("CDKN2A", "ENST00000304494", "G93W"));
        FEATURE_BLACKLIST.add(oncoKb("CDKN2A", "ENST00000304494", "V118D"));
        FEATURE_BLACKLIST.add(oncoKb("FLT3", "ENST00000241453", "L611_E612insCSSDNEYFYVDFREYEYDLKWEFPRENL"));
        FEATURE_BLACKLIST.add(oncoKb("FLT3", "ENST00000241453", "E612_F613insGYVDFREYEYDLKWEFRPRENLEF"));
        FEATURE_BLACKLIST.add(oncoKb("IGF1R", "ENST00000268035", "G119T"));
        FEATURE_BLACKLIST.add(oncoKb("IGF1R", "ENST00000268035", "G1125A"));
        FEATURE_BLACKLIST.add(oncoKb("IGF1R", "ENST00000268035", "A1374V"));
        FEATURE_BLACKLIST.add(oncoKb("JAK1", "ENST00000342505", "G1079D"));
        FEATURE_BLACKLIST.add(oncoKb("JAK1", "ENST00000342505", "G871E"));
        FEATURE_BLACKLIST.add(oncoKb("KIT", "ENST00000288135", "A504_Y505ins"));
        FEATURE_BLACKLIST.add(oncoKb("KIT", "ENST00000288135", "D814V"));
        FEATURE_BLACKLIST.add(oncoKb("MAP3K1", "ENST00000399503", "T1481fs"));
        FEATURE_BLACKLIST.add(oncoKb("PMS2", "ENST00000265849", "E541K"));
        FEATURE_BLACKLIST.add(oncoKb("POLE", "ENST00000320574", "S279Y"));
        //  - This one comes from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4155059/ and probably should be L1237F.
        FEATURE_BLACKLIST.add(oncoKb("RAD50", "ENST00000265335", "L1273F"));
        FEATURE_BLACKLIST.add(oncoKb("RAD51D", "ENST00000335858", "S257delinsK"));
        FEATURE_BLACKLIST.add(oncoKb("RXRA", "ENST00000481739", "S247F"));
        FEATURE_BLACKLIST.add(oncoKb("RXRA", "ENST00000481739", "S247Y"));
        FEATURE_BLACKLIST.add(oncoKb("SMAD4", "ENST00000342988", "D357Y"));
        FEATURE_BLACKLIST.add(oncoKb("SOX9", "ENST00000245479", "F12L"));
        FEATURE_BLACKLIST.add(oncoKb("SOX9", "ENST00000245479", "A19V"));
        FEATURE_BLACKLIST.add(oncoKb("SOX9", "ENST00000245479", "H65Y"));
        FEATURE_BLACKLIST.add(oncoKb("TERT", "ENST00000310581", "C228T"));
        FEATURE_BLACKLIST.add(oncoKb("TERT", "ENST00000310581", "C250T"));
        FEATURE_BLACKLIST.add(oncoKb("TMPRSS2", "ENST00000398585", "M160V"));

        // Variants that probably have the wrong transcripts configured by OncoKB
        FEATURE_BLACKLIST.add(oncoKb("CASP8", "ENST00000358485", "G325A"));
        FEATURE_BLACKLIST.add(oncoKb("EZH2", "ENST00000320356", "A677G"));
        FEATURE_BLACKLIST.add(oncoKb("FBXW7", "ENST00000281708", "R482Q"));
        FEATURE_BLACKLIST.add(oncoKb("FGFR2", "ENST00000358487", "K525E"));
        FEATURE_BLACKLIST.add(oncoKb("FGFR2", "ENST00000358487", "T730S"));
        FEATURE_BLACKLIST.add(oncoKb("FGFR2", "ENST00000358487", "V755I"));
        FEATURE_BLACKLIST.add(oncoKb("GNAS", "ENST00000371085", "R844H"));
        FEATURE_BLACKLIST.add(oncoKb("JAK2", "ENST00000381652", "R277K"));
        FEATURE_BLACKLIST.add(oncoKb("KRAS", "ENST00000256078", "D153V"));
        FEATURE_BLACKLIST.add(oncoKb("NF1", "ENST00000358273", "R1391G"));
        FEATURE_BLACKLIST.add(oncoKb("NF1", "ENST00000358273", "R1391S"));
        FEATURE_BLACKLIST.add(oncoKb("NF1", "ENST00000358273", "V1398D"));
        FEATURE_BLACKLIST.add(oncoKb("NF1", "ENST00000358273", "K1423E"));
        FEATURE_BLACKLIST.add(oncoKb("NF1", "ENST00000358273", "K1436Q"));
        FEATURE_BLACKLIST.add(oncoKb("NF1", "ENST00000358273", "S1463F"));
        FEATURE_BLACKLIST.add(oncoKb("NKX2-1", "ENST00000354822", "A339V"));
        FEATURE_BLACKLIST.add(oncoKb("PRDM1", "ENST00000369096", "P48R"));
        FEATURE_BLACKLIST.add(oncoKb("PRDM1", "ENST00000369096", "P48T"));
        FEATURE_BLACKLIST.add(oncoKb("PRDM1", "ENST00000369096", "Y149D"));
        FEATURE_BLACKLIST.add(oncoKb("PRDM1", "ENST00000369096", "C569Y"));
        FEATURE_BLACKLIST.add(oncoKb("PTPRT", "ENST00000373198", "T844M"));
        FEATURE_BLACKLIST.add(oncoKb("PTPRT", "ENST00000373198", "D927G"));
        FEATURE_BLACKLIST.add(oncoKb("PTPRT", "ENST00000373198", "V995M"));
        FEATURE_BLACKLIST.add(oncoKb("PTPRT", "ENST00000373198", "A1022E"));
        FEATURE_BLACKLIST.add(oncoKb("PTPRT", "ENST00000373198", "R1040L"));
        FEATURE_BLACKLIST.add(oncoKb("RAC1", "ENST00000356142", "C157Y"));
        FEATURE_BLACKLIST.add(oncoKb("RBM10", "ENST00000329236", "V354E"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "G42R"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "H78Q"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "R80C"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "K83E"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "K83N"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "Y113*"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "A122*"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "R139G"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "D171G"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "D171N"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "P173S"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "R174*"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "R174Q"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "R177*"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "R177Q"));
        FEATURE_BLACKLIST.add(oncoKb("TGFBR2", "ENST00000359013", "V419L"));
        FEATURE_BLACKLIST.add(oncoKb("TGFBR2", "ENST00000359013", "P525L"));
        FEATURE_BLACKLIST.add(oncoKb("TGFBR2", "ENST00000359013", "E526Q"));
        FEATURE_BLACKLIST.add(oncoKb("TGFBR2", "ENST00000359013", "R537P"));
        FEATURE_BLACKLIST.add(oncoKb("TMPRSS2", "ENST00000398585", "T75M"));

        // Variants that are simply invalid protein annotations
        FEATURE_BLACKLIST.add(oncoKb("BRAF", "ENST00000288602", "V600D_K601insFGLAT"));
        FEATURE_BLACKLIST.add(oncoKb("CARD11", "ENST00000396946", "L225LI"));
        FEATURE_BLACKLIST.add(oncoKb("RIT1", "ENST00000368323", "TA83del"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "S291fsX300"));
        FEATURE_BLACKLIST.add(oncoKb("RUNX1", "ENST00000300305", "S70fsX93"));
    }

    @NotNull
    private static FeatureCurationKey civic(@NotNull String gene, @Nullable String transcript, @NotNull String featureName) {
        return new FeatureCurationKey(ViccSource.CIVIC, gene, transcript, featureName);
    }

    @NotNull
    private static FeatureCurationKey cgi(@NotNull String gene, @Nullable String transcript, @NotNull String featureName) {
        return new FeatureCurationKey(ViccSource.CGI, gene, transcript, featureName);
    }

    @NotNull
    private static FeatureCurationKey jax(@NotNull String gene, @Nullable String transcript, @NotNull String featureName) {
        return new FeatureCurationKey(ViccSource.JAX, gene, transcript, featureName);
    }

    @NotNull
    private static FeatureCurationKey oncoKb(@NotNull String gene, @Nullable String transcript, @NotNull String featureName) {
        return new FeatureCurationKey(ViccSource.ONCOKB, gene, transcript, featureName);
    }

    @NotNull
    private static FeatureCurationValues curation(@NotNull String geneSymbol, @NotNull String featureName) {
        return ImmutableFeatureCurationValues.builder().geneSymbol(geneSymbol).featureName(featureName).build();
    }
}
