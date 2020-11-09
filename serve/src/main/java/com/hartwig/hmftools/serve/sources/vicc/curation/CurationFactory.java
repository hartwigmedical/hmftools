package com.hartwig.hmftools.serve.sources.vicc.curation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class CurationFactory {

    static final Map<CurationKey, String> FEATURE_NAME_MAPPINGS = Maps.newHashMap();

    static final Set<CurationKey> FEATURE_BLACKLIST = Sets.newHashSet();

    static {
        populateCGICuration();
        populateCIViCCuration();
        populateJaxCuration();
        populateOncoKBCuration();
    }

    private static void populateCGICuration() {
        // These are not wrong but lead to inconsistencies downstream
        FEATURE_NAME_MAPPINGS.put(cgiKey("PIK3R1", null, "PIK3R1 p.E439delE"), "PIK3R1 p.E439del");
        FEATURE_NAME_MAPPINGS.put(cgiKey("PIK3R1", null, "PIK3R1 p.D560_S565delDKRMNS"), "PIK3R1 p.D560_S565del");
        FEATURE_NAME_MAPPINGS.put(cgiKey("PIK3R1", null, "PIK3R1 p.T576delT"), "PIK3R1 p.T576del");

        // These fusions need to be swapped around.
        FEATURE_NAME_MAPPINGS.put(cgiKey("ABL1", null, "ABL1-BCR fusion"), "BCR-ABL1 fusion");
        FEATURE_NAME_MAPPINGS.put(cgiKey("BCR", null, "ABL1-BCR fusion"), "BCR-ABL1 fusion");
        FEATURE_NAME_MAPPINGS.put(cgiKey("COL1A1", null, "PDGFB-COL1A1 fusion"), "COL1A1-PDGFB fusion");
        FEATURE_NAME_MAPPINGS.put(cgiKey("PDGFB", null, "PDGFB-COL1A1 fusion"), "COL1A1-PDGFB fusion");
        FEATURE_NAME_MAPPINGS.put(cgiKey("FIP1L1", null, "PDGFRA-FIP1L1 fusion"), "FIP1L1-PDGFRA fusion");
        FEATURE_NAME_MAPPINGS.put(cgiKey("PDGFRA", null, "PDGFRA-FIP1L1 fusion"), "FIP1L1-PDGFRA fusion");
        FEATURE_NAME_MAPPINGS.put(cgiKey("PRKCH", null, "PRKCH amplification + ABL1-BCR fusion"), "PRKCH amplification + BCR-ABL1 fusion");

        // These fusions are on synonym genes
        FEATURE_NAME_MAPPINGS.put(cgiKey("BRD4", null, "BRD4-C15orf55 fusion"), "BRD4-NUTM1 fusion");
        FEATURE_NAME_MAPPINGS.put(cgiKey("C15orf55", null, "BRD4-C15orf55 fusion"), "BRD4-NUTM1 fusion");

        // Improve for consistency
        FEATURE_NAME_MAPPINGS.put(cgiKey("MET", null, "MET (Y1230C;Y1235D)"), "MET (Y1230C,Y1235D)");
    }

    private static void populateCIViCCuration() {
        // These are not wrong but transvar can't resolve generic start lost variants so we update to something that transvar can interpret.
        FEATURE_NAME_MAPPINGS.put(civicKey("VHL", "ENST00000256474", "M1? (c.3G>A)"), "M1I (c.3G>A)");
        FEATURE_NAME_MAPPINGS.put(civicKey("VHL", "ENST00000256474", "M1? (c.1-1_20del21)"), "M1I (c.1-1_20del21)");

        // Protein annotation is not formally correct, update to make it correct.
        FEATURE_NAME_MAPPINGS.put(civicKey("EGFR", "ENST00000275493", "V769_770insASV"), "V769_D770insASV");
        FEATURE_NAME_MAPPINGS.put(civicKey("ERBB2", "ENST00000269571", "M774INSAYVM"), "M774_A775INSAYVM");

        // Fusion needs to be flipped around
        FEATURE_NAME_MAPPINGS.put(civicKey("BRAF", "ENST00000288602", "BRAF-CUL1"), "CUL1-BRAF");

        // Fusions where gene is a synonym of our gene
        FEATURE_NAME_MAPPINGS.put(civicKey("FGFR1", null, "ZNF198-FGFR1"), "ZMYM2-FGFR1");
        FEATURE_NAME_MAPPINGS.put(civicKey("KMT2A", null, "MLL-MLLT3"), "KMT2A-MLLT3");

        // Fusions are not correct
        FEATURE_NAME_MAPPINGS.put(civicKey("ABL1", "ENST00000318560", "BCR-ABL G398R"), "BCR-ABL1 G398R");
        FEATURE_NAME_MAPPINGS.put(civicKey("ABL1", "ENST00000372348", "BCR-ABL T334I"), "BCR-ABL1 T334I");
        FEATURE_NAME_MAPPINGS.put(civicKey("ABL1", "ENST00000318560", "BCR-ABL E255K"), "BCR-ABL1 E255K");
        FEATURE_NAME_MAPPINGS.put(civicKey("ABL1", "ENST00000318560", "BCR-ABL F317L"), "BCR-ABL1 F317L");
        FEATURE_NAME_MAPPINGS.put(civicKey("ABL1", "ENST00000318560", "BCR-ABL F486S"), "BCR-ABL1 F486S");
        FEATURE_NAME_MAPPINGS.put(civicKey("ABL1", "ENST00000305877", "BCR-ABL"), "BCR-ABL1");

        // Variants that don't exist
        FEATURE_BLACKLIST.add(civicKey("NOTCH1", "ENST00000277541", "S2275FS"));
        FEATURE_BLACKLIST.add(civicKey("NOTCH1", "ENST00000277541", "V2444FS"));
        FEATURE_BLACKLIST.add(civicKey("VHL", null, "L132fs (c.395delA)"));
        FEATURE_BLACKLIST.add(civicKey("GNAS", "ENST00000371100", "T393C"));
        FEATURE_BLACKLIST.add(civicKey("TERT", "ENST00000310581", "C228T"));
        FEATURE_BLACKLIST.add(civicKey("VHL", null, "P81delRVV (c.243_251delGCGCGTCGT)"));
        FEATURE_BLACKLIST.add(civicKey("VHL", null, "A121I (c.364_365GC>AT)"));
        FEATURE_BLACKLIST.add(civicKey("VHL", null, "N167T (c.392A>C)"));
        FEATURE_BLACKLIST.add(civicKey("VHL", null, "X214L (c.641G>T)"));

        // Variant only possible through an MNV which spans an intron, and inconsistent with its coding impact.
        FEATURE_BLACKLIST.add(civicKey("VHL", null, "V155C (c.462delA)"));

        // Variant unlikely to be real as it spans multiple exons
        FEATURE_BLACKLIST.add(civicKey("VHL", null, "G114dup (c.342dupGGT)"));

        // Variant hard to interpret as it crossing exonic boundary (in protein space)
        FEATURE_BLACKLIST.add(civicKey("VHL", null, "G114fs (c.341delCGTTTCCAACAATTTCTCGGTGT)"));

        // Variant hard to interpret in protein space as it changes an amino acid followed by a frameshift.
        FEATURE_BLACKLIST.add(civicKey("VHL", "ENST00000256474", "M54IFS (c.162_166delGGAGG)"));

        // Variant has an inconsistent formatting. Not sure how to correct
        FEATURE_BLACKLIST.add(civicKey("PDGFRA", "ENST00000257290", "DI842-843VM"));
    }

    private static void populateJaxCuration() {
        // Update protein annotation to be correct (should be capitalized).
        FEATURE_NAME_MAPPINGS.put(jaxKey("PIK3CA", null, "PIK3CA E545k "), "PIK3CA E545K ");

        // These mappings are to work around the missing transcripts in JAX.
        // We map every mutation that appears on both the canonical + non-canonical form to its canonical form in ensembl.
        FEATURE_NAME_MAPPINGS.put(jaxKey("EZH2", null, "EZH2 Y641F "), "EZH2 Y646F ");
        FEATURE_NAME_MAPPINGS.put(jaxKey("EZH2", null, "EZH2 Y641H "), "EZH2 Y646H ");
        FEATURE_NAME_MAPPINGS.put(jaxKey("EZH2", null, "EZH2 Y641N "), "EZH2 Y646N ");
        FEATURE_NAME_MAPPINGS.put(jaxKey("EZH2", null, "EZH2 Y641S "), "EZH2 Y646S ");
        FEATURE_NAME_MAPPINGS.put(jaxKey("FGFR2", null, "FGFR2 V564I "), "FGFR2 V565I ");
        FEATURE_NAME_MAPPINGS.put(jaxKey("FGFR3", null, "FGFR3 Y373C "), "FGFR3 Y375C ");
        FEATURE_NAME_MAPPINGS.put(jaxKey("FGFR3", null, "FGFR3 K650E "), "FGFR3 K652E ");
        FEATURE_NAME_MAPPINGS.put(jaxKey("MET", null, "MET L1195V "), "MET L1213V ");
        FEATURE_NAME_MAPPINGS.put(jaxKey("MET", null, "MET Y1230H "), "MET Y1248H ");
        FEATURE_NAME_MAPPINGS.put(jaxKey("MET", null, "MET M1250T "), "MET M1268T ");

        // These mappings are identical but used concurrently. Confirmed to be replaced by I740_K745dup
        FEATURE_NAME_MAPPINGS.put(jaxKey("EGFR", null, "EGFR I744_K745insKIPVAI "), "EGFR K745_E746insIPVAIK ");
        FEATURE_NAME_MAPPINGS.put(jaxKey("KIT", null, "KIT V559del "), "KIT V560del ");

        // The below variants in FLT3 are from a paper where an additional R was added in the ref sequence, shifting all AAs by one position.
        // This has been corrected in current live CKB.
        FEATURE_BLACKLIST.add(jaxKey("FLT3", null, "FLT3 L611_E612insCSSDNEYFYVDFREYEYDLKWEFPRENL "));
        FEATURE_BLACKLIST.add(jaxKey("FLT3", null, "FLT3 E612_F613insGYVDFREYEYDLKWEFRPRENLEF "));

        // The transcript that should have this mutation (ENST00000507379) is annotated as 3' truncated with only 1135 AAs in ensembl)
        FEATURE_BLACKLIST.add(jaxKey("APC", null, "APC S1197* "));

        // The below is pending investigation by JAX, possibly a mistake by the paper.
        FEATURE_BLACKLIST.add(jaxKey("PIK3CA", null, "PIK3CA R425L "));

        // Variant that doesn't seem to exist.
        FEATURE_BLACKLIST.add(jaxKey("PTEN", null, "PTEN Y86fs "));

        // Variants spanning multiple exons in our interpretation
        FEATURE_BLACKLIST.add(jaxKey("PTEN", null, "PTEN R55fs*1 "));
        FEATURE_BLACKLIST.add(jaxKey("PTEN", null, "PTEN Y27fs*1 "));

        // Variant hard to interpret as it crossing exonic boundary
        FEATURE_BLACKLIST.add(jaxKey("VHL", null, "VHL V155fs "));
    }

    private static void populateOncoKBCuration() {
        // The dashes in these keys are not the usual dashes!
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("FGFR1", "ENST00000425967", "ERLIN2–FGFR1 Fusion"), "ERLIN2-FGFR1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("FGFR2", "ENST00000358487", "FGFR2–PPHLN1 Fusion"), "FGFR2-PPHLN1 Fusion");

        // The spaces are not consistent
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("FGFR3", "ENST00000260795", "FGFR3 - BAIAP2L1 Fusion"), "FGFR3-BAIAP2L1 Fusion");

        // These fusions need to be flipped.
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("CCND1", "ENST00000227507", "CCND1-IGH Fusion"), "IGH-CCND1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("ROS1", "ENST00000368508", "ROS1-CD74 Fusion"), "CD74-ROS1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("NSD1", "ENST00000439151", "NSD1-NUP98 Fusion"), "NUP98-NSD1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("RET", "ENST00000355710", "RET-CCDC6 Fusion"), "CCDC6-RET Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("EP300", "ENST00000263253", "EP300-MOZ Fusion"), "KAT6A-EP300 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("EP300", "ENST00000263253", "EP300-MLL Fusion"), "KMT2A-EP300 Fusion");

        // These fusions are incorrect but recoverable
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("PAX8", "ENST00000263334", "PAX8-PPARγ Fusion"), "PAX8-PPARG Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("ALK", "ENST00000389048", "NPM-ALK Fusion"), "NPM1-ALK Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("BRD4", "ENST00000263377", "BRD4-NUT Fusion"), "BRD4-NUTM1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("NOTCH1", "ENST00000277541", "SEC16A1-NOTCH1 Fusion"), "SEC16A-NOTCH1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("FGFR1", "ENST00000425967", "FGFR1OP1-FGFR1 Fusion"), "FGFR1OP-FGFR1 Fusion");

        // These fusions are on synonym genes
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("FGFR1", "ENST00000425967", "ZNF198-FGFR1 Fusion"), "ZMYM2-FGFR1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("TET1", "ENST00000373644", "MLL-TET1 Fusion"), "KMT2A-TET1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("ROS1", "ENST00000368508", "FIG-ROS1 Fusion"), "GOPC-ROS1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("PDGFRB", "ENST00000261799", "GPIAP1-PDGFRB Fusion"), "CAPRIN1-PDGFRB Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("PDGFRB", "ENST00000261799", "KIAA1509-PDGFRB Fusion"), "CCDC88C-PDGFRB Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("JAK2", "ENST00000381652", "TEL-JAK2 Fusion"), "ETV6-JAK2 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("RUNX1", "ENST00000300305", "TEL-RUNX1 Fusion"), "ETV6-RUNX1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("RUNX1", "ENST00000300305", "RUNX1-EVI1 Fusion"), "RUNX1-MECOM Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("FGFR1", "ENST00000425967", "CEP110-FGFR1 Fusion"), "CNTRL-FGFR1 Fusion");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("FGFR2", "ENST00000358487", "FGFR2-KIAA1967 Fusion"), "FGFR2-CCAR2 Fusion");

        // These are inconsistent or improperly aligned variants.
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("BRAF", "ENST00000288602", "T599insTT"), "T599_V600insTT");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("EGFR", "ENST00000275493", "E746_T751insIP"), "E746_L747insIP");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("EGFR", "ENST00000275493", "H773insLGNP"), "H773_V774insLGNP");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("EPAS1", "ENST00000263734", "533_534del"), "I533_P534del");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("EPAS1", "ENST00000263734", "534_536del"), "P534_D536del");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("KIT", "ENST00000288135", "V559del"), "V560del");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("PTEN", "ENST00000371953", "I32del"), "I33del");
        FEATURE_NAME_MAPPINGS.put(oncoKbKey("RIT1", "ENST00000368323", "T76insTLDT"), "T76_A77insTLDT");

        // Fusions that we don't know what gene they are on
        FEATURE_BLACKLIST.add(oncoKbKey("NKX2-1", "ENST00000354822", "TRA-NKX2-1 Fusion"));
        FEATURE_BLACKLIST.add(oncoKbKey("NTRK1", "ENST00000524377", "Delta-NTRK1 Fusion"));

        // Variants are unlikely as they span multiple exons (and hence are more fusions than inframes)
        FEATURE_BLACKLIST.add(oncoKbKey("PDGFRA", "ENST00000257290", "E311_K312del"));
        FEATURE_BLACKLIST.add(oncoKbKey("ETV6", "ENST00000396373", "385_418del"));

        // Variants that we couldn't quite figure out
        FEATURE_BLACKLIST.add(oncoKbKey("KIT", "ENST00000288135", "T574insTQLPYD"));

        // Variants that don't exist
        //  - Below is called a "silent promoter" according to https://pubmed.ncbi.nlm.nih.gov/11606402/
        FEATURE_BLACKLIST.add(oncoKbKey("APC", "ENST00000257430", "A290T"));
        FEATURE_BLACKLIST.add(oncoKbKey("AR", "ENST00000374690", "F876L"));
        FEATURE_BLACKLIST.add(oncoKbKey("BRCA1", "ENST00000357654", "T47D"));
        FEATURE_BLACKLIST.add(oncoKbKey("BRCA1", "ENST00000357654", "E60L"));
        FEATURE_BLACKLIST.add(oncoKbKey("BRCA2", "ENST00000380152", "S1670A"));
        FEATURE_BLACKLIST.add(oncoKbKey("BRCA2", "ENST00000380152", "R2304C"));
        FEATURE_BLACKLIST.add(oncoKbKey("BRCA2", "ENST00000380152", "N2829R"));
        FEATURE_BLACKLIST.add(oncoKbKey("CARD11", "ENST00000396946", "G116S"));
        FEATURE_BLACKLIST.add(oncoKbKey("CARD11", "ENST00000396946", "F123I"));
        FEATURE_BLACKLIST.add(oncoKbKey("CARD11", "ENST00000396946", "E127G"));
        FEATURE_BLACKLIST.add(oncoKbKey("CARM1", "ENST00000327064", "S217A"));
        FEATURE_BLACKLIST.add(oncoKbKey("CARM1", "ENST00000327064", "S217C"));
        FEATURE_BLACKLIST.add(oncoKbKey("CARM1", "ENST00000327064", "S217E"));
        FEATURE_BLACKLIST.add(oncoKbKey("CARM1", "ENST00000327064", "S217T"));
        FEATURE_BLACKLIST.add(oncoKbKey("CASP8", "ENST00000358485", "C248T"));
        FEATURE_BLACKLIST.add(oncoKbKey("CCND3", "ENST00000372991", "T286A"));
        FEATURE_BLACKLIST.add(oncoKbKey("CDK12", "ENST00000447079", "K765R"));
        FEATURE_BLACKLIST.add(oncoKbKey("CDK12", "ENST00000447079", "D887N"));
        FEATURE_BLACKLIST.add(oncoKbKey("CDKN2A", "ENST00000304494", "P73L"));
        FEATURE_BLACKLIST.add(oncoKbKey("CDKN2A", "ENST00000304494", "R79P"));
        FEATURE_BLACKLIST.add(oncoKbKey("CDKN2A", "ENST00000304494", "G93W"));
        FEATURE_BLACKLIST.add(oncoKbKey("CDKN2A", "ENST00000304494", "V118D"));
        FEATURE_BLACKLIST.add(oncoKbKey("FLT3", "ENST00000241453", "L611_E612insCSSDNEYFYVDFREYEYDLKWEFPRENL"));
        FEATURE_BLACKLIST.add(oncoKbKey("FLT3", "ENST00000241453", "E612_F613insGYVDFREYEYDLKWEFRPRENLEF"));
        FEATURE_BLACKLIST.add(oncoKbKey("IGF1R", "ENST00000268035", "G119T"));
        FEATURE_BLACKLIST.add(oncoKbKey("IGF1R", "ENST00000268035", "G1125A"));
        FEATURE_BLACKLIST.add(oncoKbKey("IGF1R", "ENST00000268035", "A1374V"));
        FEATURE_BLACKLIST.add(oncoKbKey("JAK1", "ENST00000342505", "G1079D"));
        FEATURE_BLACKLIST.add(oncoKbKey("JAK1", "ENST00000342505", "G871E"));
        FEATURE_BLACKLIST.add(oncoKbKey("KIT", "ENST00000288135", "A504_Y505ins"));
        FEATURE_BLACKLIST.add(oncoKbKey("KIT", "ENST00000288135", "D814V"));
        FEATURE_BLACKLIST.add(oncoKbKey("MAP3K1", "ENST00000399503", "T1481fs"));
        FEATURE_BLACKLIST.add(oncoKbKey("PMS2", "ENST00000265849", "E541K"));
        FEATURE_BLACKLIST.add(oncoKbKey("POLE", "ENST00000320574", "S279Y"));
        //  - This one comes from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4155059/ and probably should be L1237F.
        FEATURE_BLACKLIST.add(oncoKbKey("RAD50", "ENST00000265335", "L1273F"));
        FEATURE_BLACKLIST.add(oncoKbKey("RAD51D", "ENST00000335858", "S257delinsK"));
        FEATURE_BLACKLIST.add(oncoKbKey("RXRA", "ENST00000481739", "S247F"));
        FEATURE_BLACKLIST.add(oncoKbKey("RXRA", "ENST00000481739", "S247Y"));
        FEATURE_BLACKLIST.add(oncoKbKey("SMAD4", "ENST00000342988", "D357Y"));
        FEATURE_BLACKLIST.add(oncoKbKey("SOX9", "ENST00000245479", "F12L"));
        FEATURE_BLACKLIST.add(oncoKbKey("SOX9", "ENST00000245479", "A19V"));
        FEATURE_BLACKLIST.add(oncoKbKey("SOX9", "ENST00000245479", "H65Y"));
        FEATURE_BLACKLIST.add(oncoKbKey("TERT", "ENST00000310581", "C228T"));
        FEATURE_BLACKLIST.add(oncoKbKey("TERT", "ENST00000310581", "C250T"));
        FEATURE_BLACKLIST.add(oncoKbKey("TMPRSS2", "ENST00000398585", "M160V"));

        // Variants that probably have the wrong transcripts configured by OncoKB
        FEATURE_BLACKLIST.add(oncoKbKey("CASP8", "ENST00000358485", "G325A"));
        FEATURE_BLACKLIST.add(oncoKbKey("EZH2", "ENST00000320356", "A677G"));
        FEATURE_BLACKLIST.add(oncoKbKey("FBXW7", "ENST00000281708", "R482Q"));
        FEATURE_BLACKLIST.add(oncoKbKey("FGFR2", "ENST00000358487", "K525E"));
        FEATURE_BLACKLIST.add(oncoKbKey("FGFR2", "ENST00000358487", "T730S"));
        FEATURE_BLACKLIST.add(oncoKbKey("FGFR2", "ENST00000358487", "V755I"));
        FEATURE_BLACKLIST.add(oncoKbKey("GNAS", "ENST00000371085", "R844H"));
        FEATURE_BLACKLIST.add(oncoKbKey("JAK2", "ENST00000381652", "R277K"));
        FEATURE_BLACKLIST.add(oncoKbKey("KRAS", "ENST00000256078", "D153V"));
        FEATURE_BLACKLIST.add(oncoKbKey("NF1", "ENST00000358273", "R1391G"));
        FEATURE_BLACKLIST.add(oncoKbKey("NF1", "ENST00000358273", "R1391S"));
        FEATURE_BLACKLIST.add(oncoKbKey("NF1", "ENST00000358273", "V1398D"));
        FEATURE_BLACKLIST.add(oncoKbKey("NF1", "ENST00000358273", "K1423E"));
        FEATURE_BLACKLIST.add(oncoKbKey("NF1", "ENST00000358273", "K1436Q"));
        FEATURE_BLACKLIST.add(oncoKbKey("NF1", "ENST00000358273", "S1463F"));
        FEATURE_BLACKLIST.add(oncoKbKey("NKX2-1", "ENST00000354822", "A339V"));
        FEATURE_BLACKLIST.add(oncoKbKey("PRDM1", "ENST00000369096", "P48R"));
        FEATURE_BLACKLIST.add(oncoKbKey("PRDM1", "ENST00000369096", "P48T"));
        FEATURE_BLACKLIST.add(oncoKbKey("PRDM1", "ENST00000369096", "Y149D"));
        FEATURE_BLACKLIST.add(oncoKbKey("PRDM1", "ENST00000369096", "C569Y"));
        FEATURE_BLACKLIST.add(oncoKbKey("PTPRT", "ENST00000373198", "T844M"));
        FEATURE_BLACKLIST.add(oncoKbKey("PTPRT", "ENST00000373198", "D927G"));
        FEATURE_BLACKLIST.add(oncoKbKey("PTPRT", "ENST00000373198", "V995M"));
        FEATURE_BLACKLIST.add(oncoKbKey("PTPRT", "ENST00000373198", "A1022E"));
        FEATURE_BLACKLIST.add(oncoKbKey("PTPRT", "ENST00000373198", "R1040L"));
        FEATURE_BLACKLIST.add(oncoKbKey("RAC1", "ENST00000356142", "C157Y"));
        FEATURE_BLACKLIST.add(oncoKbKey("RBM10", "ENST00000329236", "V354E"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "G42R"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "H78Q"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "R80C"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "K83E"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "K83N"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "Y113*"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "A122*"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "R139G"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "D171G"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "D171N"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "P173S"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "R174*"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "R174Q"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "R177*"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "R177Q"));
        FEATURE_BLACKLIST.add(oncoKbKey("TGFBR2", "ENST00000359013", "V419L"));
        FEATURE_BLACKLIST.add(oncoKbKey("TGFBR2", "ENST00000359013", "P525L"));
        FEATURE_BLACKLIST.add(oncoKbKey("TGFBR2", "ENST00000359013", "E526Q"));
        FEATURE_BLACKLIST.add(oncoKbKey("TGFBR2", "ENST00000359013", "R537P"));
        FEATURE_BLACKLIST.add(oncoKbKey("TMPRSS2", "ENST00000398585", "T75M"));

        // Variants that are simply invalid protein annotations
        FEATURE_BLACKLIST.add(oncoKbKey("BRAF", "ENST00000288602", "V600D_K601insFGLAT"));
        FEATURE_BLACKLIST.add(oncoKbKey("CARD11", "ENST00000396946", "L225LI"));
        FEATURE_BLACKLIST.add(oncoKbKey("RIT1", "ENST00000368323", "TA83del"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "S291fsX300"));
        FEATURE_BLACKLIST.add(oncoKbKey("RUNX1", "ENST00000300305", "S70fsX93"));
    }

    @NotNull
    private static CurationKey civicKey(@NotNull String gene, @Nullable String transcript, @NotNull String featureName) {
        return new CurationKey(ViccSource.CIVIC, gene, transcript, featureName);
    }

    @NotNull
    private static CurationKey cgiKey(@NotNull String gene, @Nullable String transcript, @NotNull String featureName) {
        return new CurationKey(ViccSource.CGI, gene, transcript, featureName);
    }

    @NotNull
    private static CurationKey jaxKey(@NotNull String gene, @Nullable String transcript, @NotNull String featureName) {
        return new CurationKey(ViccSource.JAX, gene, transcript, featureName);
    }

    @NotNull
    private static CurationKey oncoKbKey(@NotNull String gene, @Nullable String transcript, @NotNull String featureName) {
        return new CurationKey(ViccSource.ONCOKB, gene, transcript, featureName);
    }

    private CurationFactory() {
    }
}
