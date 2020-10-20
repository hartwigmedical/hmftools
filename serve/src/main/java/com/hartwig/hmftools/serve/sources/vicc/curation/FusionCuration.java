package com.hartwig.hmftools.serve.sources.vicc.curation;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class FusionCuration {

    public FusionCuration(){
    }

    @NotNull
    public static String curatedFusions(@NotNull String fusion) {
        if (fusion.equals("ZNF198-FGFR1")) {
            fusion = "ZMYM2-FGFR1";
        } else if (fusion.equals("NPM-ALK")) {
            fusion = "NPM1-ALK";
        } else if (fusion.equals("ABL1-BCR")) {
            fusion = "BCR-ABL1";
        } else if (fusion.equals("ROS1-CD74")) {
            fusion = "CD74-ROS1";
        } else if (fusion.equals("RET-CCDC6")) {
            fusion = "CCDC6-RET";
        } else if (fusion.equals("EP300-MOZ")) {
            fusion = "KAT6A-EP300";
        } else if (fusion.equals("EP300-MLL")) {
            fusion = "KMT2A-EP300";
        } else if (fusion.equals("BRD4-NUT")) {
            fusion = "BRD4-NUTM1";
        } else if (fusion.equals("CEP110-FGFR1")) {
            fusion = "CNTRL-FGFR1";
        } else if (fusion.equals("FGFR2-KIAA1967")) {
            fusion = "FGFR2-CCAR2";
        } else if (fusion.equals("FIG-ROS1")) {
            fusion = "GOPC-ROS1";
        } else if (fusion.equals("GPIAP1-PDGFRB")) {
            fusion = "CAPRIN1-PDGFRB";
        } else if (fusion.equals("IGL-MYC")) {
            fusion = "IGLC6-MYC";
        } else if (fusion.equals("KIAA1509-PDGFRB")) {
            fusion = "CCDC88C-PDGFRB";
        } else if (fusion.equals("MLL-TET1")) {
            fusion = "KMT2A-TET1";
        } else if (fusion.equals("PAX8-PPAR?")) {
            fusion = "PAX8-PPARA";
        } else if (fusion.equals("SEC16A1-NOTCH1")) {
            fusion = "SEC16A-NOTCH1";
        } else if (fusion.equals("TEL-JAK2")) {
            fusion = "ETV6-JAK2";
        } else if (fusion.equals("TRA-NKX2-1")) {
            fusion = "TRAC-NKX2-1";
        } else if (fusion.equals("FGFR1OP1-FGFR1")) {
            fusion = "FGFR1OP-FGFR1";
        } else if (fusion.equals("PDGFRA-FIP1L1")) {
            fusion = "FIP1L1-PDGFRA";
        } else if (fusion.equals("PDGFB-COL1A1")) {
            fusion = "COL1A1-PDGFB";
        } else if (fusion.equals("BRD4-C15orf55")) {
            fusion = "BRD4-NUTM1";
        } else if (fusion.equals("MLL-MLLT3")) {
            fusion = "KMT2A-MLLT3";
        } else if (fusion.equals("BCR-ABL")) {
            fusion = "BCR-ABL1";
        } else if (fusion.equals("ERLIN2?FGFR1")) {
            fusion = "ERLIN2-FGFR1";
        } else if (fusion.equals("FGFR2?PPHLN1")) {
            fusion = "FGFR2-PPHLN1";
        } else if (fusion.equals("FGFR3 - BAIAP2L1")) {
            fusion = "FGFR3-BAIAP2L1";
        } else if (fusion.equals("PAX8-PPAR?")) {
            fusion = "PAX8-PPARA";
        }
        //TODO determine if below is needed for removing fusions from knowledgebase
        else if (fusion.contains("IGH") || fusion.contains("IGK") || fusion.contains("TRB") || fusion.contains("Delta") || fusion.equals(
                "RET-TPCN1") || fusion.equals("PVT1-MYC") || fusion.equals("ESR1-CCDC170") || fusion.equals("BRAF-CUL1")) {
            //function.contains("Loss-of-function")
            fusion = Strings.EMPTY;
        }
        return fusion;
    }
}
