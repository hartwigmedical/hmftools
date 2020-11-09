package com.hartwig.hmftools.vicc.annotation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ComplexClassifier {

    private static final Map<String, Set<String>> COMPLEX_EVENTS_PER_GENE = createComplexEventMap();

    public static boolean isComplexEvent(@NotNull String featureName, @Nullable String gene) {
        Set<String> entriesForGene = COMPLEX_EVENTS_PER_GENE.get(gene);
        if (entriesForGene != null && entriesForGene.contains(featureName.trim())) {
            return true;
        } else if (featureName.trim().endsWith(".")) {
            // Not sure what a dot means!
            return true;
        } else if (featureName.split("\\*").length > 2) {
            // Some hotspots contain multiple stop codons.
            return true;
        } else if (featureName.contains("/")) {
            // Some hotspots appear as V600E/K
            return !featureName.toLowerCase().contains("exon");
        } else {
            // Some frameshifts also change the amino acid itself in the position of the frameshift.
            int fsLocation = featureName.indexOf("fs");
            return fsLocation > 1 && !isInteger(featureName.substring(fsLocation - 1, fsLocation));
        }
    }

    private static boolean isInteger(@NotNull String string) {
        try {
            Integer.parseInt(string);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    @NotNull
    private static Map<String, Set<String>> createComplexEventMap() {
        Map<String, Set<String>> complexEventsPerGene = Maps.newHashMap();

        Set<String> alkSet = Sets.newHashSet("ALK inframe insertion (1151T)");
        complexEventsPerGene.put("ALK", alkSet);

        Set<String> arSet = Sets.newHashSet("ARv567es", "SPLICE VARIANT 7", "AR-V7");
        complexEventsPerGene.put("AR", arSet);

        Set<String> brafSet = Sets.newHashSet("DEL 485-490", "L485_P490>Y");
        complexEventsPerGene.put("BRAF", brafSet);

        Set<String> brca1Set = Sets.newHashSet("Alu insertion");
        complexEventsPerGene.put("BRCA1", brca1Set);

        Set<String> brca2Set =
                Sets.newHashSet("V1790_A1996del", "T1815_Y1894del", "S1834_G1892del", "V1839_E1901del", "A1847_M1890del", "S1871_C1893del");
        complexEventsPerGene.put("BRCA2", brca2Set);

        Set<String> casp8Set = Sets.newHashSet("CASP8L");
        complexEventsPerGene.put("CASP8", casp8Set);

        Set<String> cebpaSet = Sets.newHashSet("N-TERMINAL FRAME SHIFT");
        complexEventsPerGene.put("CEBPA", cebpaSet);

        Set<String> ccnd1Set = Sets.newHashSet("256_286trunc");
        complexEventsPerGene.put("CCND1", ccnd1Set);

        Set<String> ccnd3Set = Sets.newHashSet("D286_L292trunc");
        complexEventsPerGene.put("CCND3", ccnd3Set);

        Set<String> chek2Set = Sets.newHashSet("1100DELC", "IVS2+1G>A");
        complexEventsPerGene.put("CHEK2", chek2Set);

        Set<String> dnmt3bSet = Sets.newHashSet("DNMT3B7");
        complexEventsPerGene.put("DNMT3B", dnmt3bSet);

        Set<String> dpydSet = Sets.newHashSet("DPYD splice donor variant", "DPYD*13 HOMOZYGOSITY", "DPYD*2A HOMOZYGOSITY");
        complexEventsPerGene.put("DPYD", dpydSet);

        Set<String> egfrSet = Sets.newHashSet("EGFR L698_S1037dup",
                "EGFR inframe deletion (30-336)",
                "EGFR inframe insertion (769-770)",
                "EGFR inframe deletion (6-273)",
                "T34_A289del",
                "EGFR T34_A289del",
                "G983_A1210del",
                "EGFR CTD");
        complexEventsPerGene.put("EGFR", egfrSet);

        Set<String> eif1axSet = Sets.newHashSet("A113_splice");
        complexEventsPerGene.put("EIF1AX", eif1axSet);

        Set<String> erbb2Set = Sets.newHashSet("P780INS", "DEL 755-759");
        complexEventsPerGene.put("ERBB2", erbb2Set);

        Set<String> fli1Set = Sets.newHashSet("EWSR1-FLI1 Type 1");
        complexEventsPerGene.put("FLI1", fli1Set);

        Set<String> flt3Set = Sets.newHashSet("FLT3-ITD", "ITD", "FLT3 internal tandem duplications");
        complexEventsPerGene.put("FLT3", flt3Set);

        Set<String> hlaaSet = Sets.newHashSet("596_619splice");
        complexEventsPerGene.put("HLA-A", hlaaSet);

        Set<String> kdm5cSet = Sets.newHashSet("M1_E165DEL");
        complexEventsPerGene.put("KDM5C", kdm5cSet);

        Set<String> kitSet = Sets.newHashSet("INTERNAL DUPLICATION",
                "E554_I571del",
                "V555_L576del",
                "K550_G592del",
                "KIT V560_Y578del",
                "KIT V560_L576del",
                "KIT inframe deletion (577-579)",
                "KIT inframe deletion (V560)");
        complexEventsPerGene.put("KIT", kitSet);

        Set<String> krasSet = Sets.newHashSet("LCS6-variant");
        complexEventsPerGene.put("KRAS", krasSet);

        Set<String> map2k1Set = Sets.newHashSet("MAP2K1 inframe deletion (56-60)");
        complexEventsPerGene.put("MAP2K1", map2k1Set);

        Set<String> metSet = Sets.newHashSet("963_D1010splice",
                "981_1028splice",
                "X1008_splice",
                "X963_splice",
                "X1006_splice",
                "X1009_splice",
                "X1007_splice");
        complexEventsPerGene.put("MET", metSet);

        Set<String> mlh1Set = Sets.newHashSet("C.790+1G>A");
        complexEventsPerGene.put("MLH1", mlh1Set);

        Set<String> notch2Set = Sets.newHashSet("2010_2471trunc", "1_2009trunc");
        complexEventsPerGene.put("NOTCH2", notch2Set);

        Set<String> ntrk1Set = Sets.newHashSet("TRKAIII Splice Variant");
        complexEventsPerGene.put("NTRK1", ntrk1Set);

        Set<String> pdgfraSet = Sets.newHashSet("PDGFRA inframe deletion (I843)", "C456_R481del", "Y375_K455del");
        complexEventsPerGene.put("PDGFRA", pdgfraSet);

        Set<String> pik3r1Set = Sets.newHashSet("X582_splice", "X475_splice", "X434_splice");
        complexEventsPerGene.put("PIK3R1", pik3r1Set);

        Set<String> poleSet = Sets.newHashSet("POLE (268-471)");
        complexEventsPerGene.put("POLE", poleSet);

        Set<String> ppm1dSet = Sets.newHashSet("422_605trunc");
        complexEventsPerGene.put("PPM1D", ppm1dSet);

        Set<String> ptch1Set = Sets.newHashSet("LOH");
        complexEventsPerGene.put("PTCH1", ptch1Set);

        Set<String> runx1Set = Sets.newHashSet("R135FSX177", "T148HFSX9");
        complexEventsPerGene.put("RUNX1", runx1Set);

        Set<String> tertSet = Sets.newHashSet("TERT promoters core", "Promoter Mutations");
        complexEventsPerGene.put("TERT", tertSet);

        Set<String> tgfbr1Set = Sets.newHashSet("TGFBR1*6A");
        complexEventsPerGene.put("TGFBR1", tgfbr1Set);

        Set<String> tp53Set = Sets.newHashSet("DNA binding domain deletions", "DNA binding domain insertions");
        complexEventsPerGene.put("TP53", tp53Set);

        Set<String> tpmtSet = Sets.newHashSet("TPMT splice acceptor variant");
        complexEventsPerGene.put("TPMT", tpmtSet);

        Set<String> tymsSet = Sets.newHashSet("5' TANDEM REPEAT");
        complexEventsPerGene.put("TYMS", tymsSet);

        Set<String> ugt1a1Set = Sets.newHashSet("UGT1A1*28", "UGT1A1*60");
        complexEventsPerGene.put("UGT1A1", ugt1a1Set);

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
        complexEventsPerGene.put("VHL", vhlSet);

        return complexEventsPerGene;
    }
}
