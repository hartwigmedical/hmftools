package com.hartwig.hmftools.serve.sources.vicc.curation;

import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class FusionCuration {

    private FusionCuration() {
    }

    @NotNull
    public static String curateFusion(@NotNull String fusion, @NotNull Feature feature) {
        if (fusion.equals("ZNF198-FGFR1")) {
            return "ZMYM2-FGFR1";
        } else if (fusion.equals("NPM-ALK")) {
            return "NPM1-ALK";
        } else if (fusion.equals("ABL1-BCR")) {
            return "BCR-ABL1";
        } else if (fusion.equals("ROS1-CD74")) {
            return "CD74-ROS1";
        } else if (fusion.equals("RET-CCDC6")) {
            return "CCDC6-RET";
        } else if (fusion.equals("EP300-MOZ")) {
            return "KAT6A-EP300";
        } else if (fusion.equals("EP300-MLL")) {
            return "KMT2A-EP300";
        } else if (fusion.equals("BRD4-NUT")) {
            return "BRD4-NUTM1";
        } else if (fusion.equals("CEP110-FGFR1")) {
            return "CNTRL-FGFR1";
        } else if (fusion.equals("FGFR2-KIAA1967")) {
            return "FGFR2-CCAR2";
        } else if (fusion.equals("FIG-ROS1")) {
            return "GOPC-ROS1";
        } else if (fusion.equals("GPIAP1-PDGFRB")) {
            return "CAPRIN1-PDGFRB";
        } else if (fusion.equals("IGL-MYC")) {
            return "IGLC6-MYC";
        } else if (fusion.equals("KIAA1509-PDGFRB")) {
            return "CCDC88C-PDGFRB";
        } else if (fusion.equals("MLL-TET1")) {
            return "KMT2A-TET1";
        } else if (fusion.equals("PAX8-PPAR?")) {
            return "PAX8-PPARA";
        } else if (fusion.equals("SEC16A1-NOTCH1")) {
            return "SEC16A-NOTCH1";
        } else if (fusion.equals("TEL-JAK2")) {
            return "ETV6-JAK2";
        } else if (fusion.equals("TRA-NKX2-1")) {
            return "TRAC-NKX2-1";
        } else if (fusion.equals("FGFR1OP1-FGFR1")) {
            return "FGFR1OP-FGFR1";
        } else if (fusion.equals("PDGFRA-FIP1L1")) {
            return "FIP1L1-PDGFRA";
        } else if (fusion.equals("PDGFB-COL1A1")) {
            return "COL1A1-PDGFB";
        } else if (fusion.equals("BRD4-C15orf55")) {
            return "BRD4-NUTM1";
        } else if (fusion.equals("MLL-MLLT3")) {
            return "KMT2A-MLLT3";
        } else if (fusion.equals("BCR-ABL")) {
            return "BCR-ABL1";
        } else if (fusion.equals("ERLIN2?FGFR1")) {
            return "ERLIN2-FGFR1";
        } else if (fusion.equals("FGFR2?PPHLN1")) {
            return "FGFR2-PPHLN1";
        } else if (fusion.equals("FGFR3 - BAIAP2L1")) {
            return "FGFR3-BAIAP2L1";
        } else if (fusion.equals("PAX8-PPAR?")) {
            return "PAX8-PPARA";
        } else if (feature.name().equals("ITD") && feature.geneSymbol().equals("FLT3")) {
            //removed the FLT3-ITD fusion
            return Strings.EMPTY;
        }
        //TODO determine if below is needed for removing fusions from knowledgebase
        else if (fusion.contains("IGH") || fusion.contains("IGK") || fusion.contains("TRB") || fusion.contains("Delta") || fusion.equals(
                "RET-TPCN1") || fusion.equals("PVT1-MYC") || fusion.equals("ESR1-CCDC170") || fusion.equals("BRAF-CUL1")) {
            //function.contains("Loss-of-function")
            return Strings.EMPTY;
        }

        return fusion;
    }
}
