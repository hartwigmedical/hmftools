package com.hartwig.hmftools.serve.vicc.curation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

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
        // Below is not wrong but leads to inconsistencies downstream
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.CGI, "PIK3R1", null, "PIK3R1 p.E439delE"), "PIK3R1 p.E439del");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.CGI, "PIK3R1", null, "PIK3R1 p.D560_S565delDKRMNS"), "PIK3R1 p.D560_S565del");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.CGI, "PIK3R1", null, "PIK3R1 p.T576delT"), "PIK3R1 p.T576del");
    }

    private static void populateCIViCCuration() {
        // Below is not wrong but transvar struggles with interpreting start lost variants,
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.CIVIC, "VHL", "ENST00000256474", "M1? (c.3G>A)"), "M1I (c.3G>A)");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.CIVIC, "VHL", "ENST00000256474", "M1? (c.1-1_20del21)"),
                "M1I (c.1-1_20del21)");

        // Below is not wrong but leads to inconsistencies downstream.
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.CIVIC, "EGFR", "ENST00000275493", "V769_770insASV"),
                "V769_D770insASV");

        // Variants implying stop lost. They are real but not handled yet in SERVE (TODO DEV-1475)
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "VHL", "ENST00000256474", "*214W (c.642A>G)"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "MLH1", "ENST00000231790", "*757L"));

        // Variants that don't exist
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "NOTCH1", "ENST00000277541", "S2275FS"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "NOTCH1", "ENST00000277541", "V2444FS"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "VHL", null, "L132fs (c.395delA)"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "GNAS", "ENST00000371100", "T393C"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "TERT", "ENST00000310581", "C228T"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "VHL", null, "P81delRVV (c.243_251delGCGCGTCGT)"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "VHL", null, "A121I (c.364_365GC>AT)"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "VHL", null, "N167T (c.392A>C)"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "VHL", null, "X214L (c.641G>T)"));

        // Variant only possible through an MNV which spans an intron
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "VHL", null, "V155C (c.462delA)"));

        // Synonymous variant, unlikely to have an impact
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "VHL", null, "R161R (c.481C>A)"));

        // Variant unlikely to be real as it spans multiple exons
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "VHL", null, "G114dup (c.342dupGGT)"));

        // Variant hard to interpret as it crossing exonic boundary
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.CIVIC, "VHL", null, "G114fs (c.341delCGTTTCCAACAATTTCTCGGTGT)"));
    }

    private static void populateJaxCuration() {
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "PIK3CA", null, "PIK3CA E545k "), "PIK3CA E545K ");

        // These mappings are to work around the missing transcripts in JAX.
        // We map every mutation that appears on both the canonical + non-canonical form to its canonical form in ensembl.
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "EZH2", null, "EZH2 Y641F "), "EZH2 Y646F ");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "EZH2", null, "EZH2 Y641H "), "EZH2 Y646H ");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "EZH2", null, "EZH2 Y641N "), "EZH2 Y646N ");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "EZH2", null, "EZH2 Y641S "), "EZH2 Y646S ");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "FGFR2", null, "FGFR2 V564I "), "FGFR2 V565I ");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "FGFR3", null, "FGFR3 Y373C "), "FGFR3 Y375C ");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "FGFR3", null, "FGFR3 K650E "), "FGFR3 K652E ");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "MET", null, "MET L1195V "), "MET L1213V ");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "MET", null, "MET Y1230H "), "MET Y1248H ");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "MET", null, "MET M1250T "), "MET M1268T ");

        // These mappings are identical but used concurrently. Confirmed to be replaced by I740_K745dup
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "EGFR", null, "EGFR I744_K745insKIPVAI "),
                "EGFR K745_E746insIPVAIK ");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.JAX, "KIT", null, "KIT V559del "), "KIT V560del ");

        // The below variants in FLT3 are from a paper where an additional R was added in the ref sequence, shifting all AAs by one position.
        // This has been corrected in current live CKB.
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.JAX, "FLT3", null, "FLT3 L611_E612insCSSDNEYFYVDFREYEYDLKWEFPRENL "));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.JAX, "FLT3", null, "FLT3 E612_F613insGYVDFREYEYDLKWEFRPRENLEF "));

        // The transcript that should have this mutation (ENST00000507379.1) is annotated as 3' truncated with only 1135 AAs in ensembl)
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.JAX, "APC", null, "APC S1197* "));

        // The below is pending investigation by JAX, possibly a mistake by the paper.
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.JAX, "PIK3CA", null, "PIK3CA R425L "));

        // Variant that doesn't seem to exist. TODO Report to JAX in case still exists.
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.JAX, "PTEN", null, "PTEN Y86fs "));

        // Variant hard to interpret as it crossing exonic boundary
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.JAX, "VHL", null, "VHL V155fs "));
    }

    private static void populateOncoKBCuration() {
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.ONCOKB, "EPAS1", "ENST00000263734", "533_534del"), "I533_P534del");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.ONCOKB, "EPAS1", "ENST00000263734", "534_536del"), "P534_D536del");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.ONCOKB, "KIT", "ENST00000288135", "V559del"), "V560del");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.ONCOKB, "PTEN", "ENST00000371953", "I32del"), "I33del");
        FEATURE_NAME_MAPPINGS.put(new CurationKey(ViccSource.ONCOKB, "EGFR", "ENST00000275493", "E746_T751insIP"), "E746_L747insIP");

        // Variant is fine, but interpretation currently not supported by SERVE
        // The unaligned insert lies outside of exonic range, but the left-aligned insert is fine
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "BRAF", "ENST00000288602", "R506_K507insVLR"));

        // Variants are unlikely as they span multiple exons (and hence are more fusions than inframes)
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PDGFRA", "ENST00000257290", "E311_K312del"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "ETV6", "ENST00000396373", "385_418del"));

        // Variants that don't exist
        //  - Below is called a "silent promoter" according to https://pubmed.ncbi.nlm.nih.gov/11606402/
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "APC", "ENST00000257430", "A290T"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "AR", "ENST00000374690", "F876L"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "BRCA1", "ENST00000357654", "T47D"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "BRCA1", "ENST00000357654", "E60L"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "BRCA2", "ENST00000380152", "S1670A"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "BRCA2", "ENST00000380152", "R2304C"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "BRCA2", "ENST00000380152", "N2829R"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CARD11", "ENST00000396946", "G116S"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CARD11", "ENST00000396946", "F123I"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CARD11", "ENST00000396946", "E127G"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CARM1", "ENST00000327064", "S217A"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CARM1", "ENST00000327064", "S217C"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CARM1", "ENST00000327064", "S217E"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CARM1", "ENST00000327064", "S217T"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CASP8", "ENST00000358485", "C248T"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CCND3", "ENST00000372991", "T286A"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CDK12", "ENST00000447079", "K765R"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CDK12", "ENST00000447079", "D887N"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CDKN2A", "ENST00000304494", "P73L"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CDKN2A", "ENST00000304494", "R79P"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CDKN2A", "ENST00000304494", "G93W"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CDKN2A", "ENST00000304494", "V118D"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB,
                "FLT3",
                "ENST00000241453",
                "L611_E612insCSSDNEYFYVDFREYEYDLKWEFPRENL"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "FLT3", "ENST00000241453", "E612_F613insGYVDFREYEYDLKWEFRPRENLEF"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "IGF1R", "ENST00000268035", "G119T"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "IGF1R", "ENST00000268035", "G1125A"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "IGF1R", "ENST00000268035", "A1374V"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "JAK1", "ENST00000342505", "G1079D"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "JAK1", "ENST00000342505", "G871E"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "KIT", "ENST00000288135", "A504_Y505ins"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "KIT", "ENST00000288135", "D814V"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "MAP3K1", "ENST00000399503", "T1481fs"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PMS2", "ENST00000265849", "E541K"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "POLE", "ENST00000320574", "S279Y"));
        //  - This one comes from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4155059/ and probably should be L1237F.
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RAD50", "ENST00000265335", "L1273F"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RAD51D", "ENST00000335858", "S257delinsK"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RXRA", "ENST00000481739", "S247F"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RXRA", "ENST00000481739", "S247Y"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "SMAD4", "ENST00000342988", "D357Y"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "SOX9", "ENST00000245479", "F12L"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "SOX9", "ENST00000245479", "A19V"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "SOX9", "ENST00000245479", "H65Y"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "TERT", "ENST00000310581", "C228T"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "TERT", "ENST00000310581", "C250T"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "TMPRSS2", "ENST00000398585", "M160V"));

        // The below probably have the wrong transcripts configured by OncoKB
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CASP8", "ENST00000358485", "G325A"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "EZH2", "ENST00000320356", "A677G"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "FBXW7", "ENST00000281708", "R482Q"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "FGFR2", "ENST00000358487", "K525E"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "FGFR2", "ENST00000358487", "T730S"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "FGFR2", "ENST00000358487", "V755I"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "GNAS", "ENST00000371085", "R844H"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "JAK2", "ENST00000381652", "R277K"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "KRAS", "ENST00000256078", "D153V"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "NF1", "ENST00000358273", "R1391G"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "NF1", "ENST00000358273", "R1391S"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "NF1", "ENST00000358273", "V1398D"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "NF1", "ENST00000358273", "K1423E"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "NF1", "ENST00000358273", "K1436Q"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "NF1", "ENST00000358273", "S1463F"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "NKX2-1", "ENST00000354822", "A339V"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PRDM1", "ENST00000369096", "P48R"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PRDM1", "ENST00000369096", "P48T"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PRDM1", "ENST00000369096", "Y149D"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PRDM1", "ENST00000369096", "C569Y"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PTPRT", "ENST00000373198", "T844M"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PTPRT", "ENST00000373198", "D927G"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PTPRT", "ENST00000373198", "V995M"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PTPRT", "ENST00000373198", "A1022E"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "PTPRT", "ENST00000373198", "R1040L"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RAC1", "ENST00000356142", "C157Y"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RBM10", "ENST00000329236", "V354E"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "G42R"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "H78Q"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "R80C"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "K83E"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "K83N"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "Y113*"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "A122*"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "R139G"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "D171G"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "D171N"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "P173S"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "R174*"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "R174Q"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "R177*"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "RUNX1", "ENST00000300305", "R177Q"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "TGFBR2", "ENST00000359013", "V419L"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "TGFBR2", "ENST00000359013", "P525L"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "TGFBR2", "ENST00000359013", "E526Q"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "TGFBR2", "ENST00000359013", "R537P"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "TMPRSS2", "ENST00000398585", "T75M"));

        // The below are simply invalid protein annotations
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "BRAF", "ENST00000288602", "V600D_K601insFGLAT"));
        FEATURE_BLACKLIST.add(new CurationKey(ViccSource.ONCOKB, "CARD11", "ENST00000396946", "L225LI"));
    }

    private CurationFactory() {
    }
}
