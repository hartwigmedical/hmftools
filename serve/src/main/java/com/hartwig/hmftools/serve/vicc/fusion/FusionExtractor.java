package com.hartwig.hmftools.serve.vicc.fusion;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class FusionExtractor {

    public FusionExtractor() {
    }

    public Map<Feature, FusionAnnotation> extractKnownFusions(@NotNull ViccEntry viccEntry) {
        Map<Feature, FusionAnnotation> fusionsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.FUSION_PAIR) {
                // TODO Implement geneUp, exonUp, geneDown, exonDown
                fusionsPerFeature.put(feature,
                        ImmutableFusionAnnotation.builder().geneUp(feature.name()).geneDown(feature.name()).build());
            }
        }
        return fusionsPerFeature;
    }

    // TODO Move PROMISCUOUS fusions to Gene level event extractor.
//    @NotNull
//    public static KnownFusions determinePromiscuousFusions(@NotNull ViccSource source, @NotNull String typeEvent, @NotNull String gene,
//            @NotNull String function) {
//        if (function.equals("Likely Loss-of-function")) {
//            gene = Strings.EMPTY;
//            typeEvent = Strings.EMPTY;
//        }
//
//        return ImmutableKnownFusions.builder()
//                .gene(gene)
//                .eventType(typeEvent)
//                .source(source.toString())
//                .sourceLink(source.toString())
//                .build();
//    }
//
//    @NotNull
//    public static KnownFusions determineKnownFusionsPairs(@NotNull ViccSource source, @NotNull String typeEvent, @NotNull String gene,
//            @NotNull String function) {
//        if (gene.equals("ZNF198-FGFR1")) {
//            gene = "ZMYM2-FGFR1";
//        } else if (gene.equals("NPM-ALK")) {
//            gene = "NPM1-ALK";
//        } else if (gene.equals("ABL1-BCR")) {
//            gene = "BCR-ABL1";
//        } else if (gene.equals("ROS1-CD74")) {
//            gene = "CD74-ROS1";
//        } else if (gene.equals("RET-CCDC6")) {
//            gene = "CCDC6-RET";
//        } else if (gene.equals("EP300-MOZ")) {
//            gene = "KAT6A-EP300";
//        } else if (gene.equals("EP300-MLL")) {
//            gene = "KMT2A-EP300";
//        } else if (gene.equals("BRD4-NUT")) {
//            gene = "BRD4-NUTM1";
//        } else if (gene.equals("CEP110-FGFR1")) {
//            gene = "CNTRL-FGFR1";
//        } else if (gene.equals("FGFR2-KIAA1967")) {
//            gene = "FGFR2-CCAR2";
//        } else if (gene.equals("FIG-ROS1")) {
//            gene = "GOPC-ROS1";
//        } else if (gene.equals("GPIAP1-PDGFRB")) {
//            gene = "CAPRIN1-PDGFRB";
//        } else if (gene.equals("IGL-MYC")) {
//            gene = "IGLC6-MYC";
//        } else if (gene.equals("KIAA1509-PDGFRB")) {
//            gene = "CCDC88C-PDGFRB";
//        } else if (gene.equals("MLL-TET1")) {
//            gene = "KMT2A-TET1";
//        } else if (gene.equals("PAX8-PPAR?")) {
//            gene = "PAX8-PPARA";
//        } else if (gene.equals("SEC16A1-NOTCH1")) {
//            gene = "SEC16A-NOTCH1";
//        } else if (gene.equals("TEL-JAK2")) {
//            gene = "ETV6-JAK2";
//        } else if (gene.equals("TRA-NKX2-1")) {
//            gene = "TRAC-NKX2-1";
//        } else if (gene.equals("FGFR1OP1-FGFR1")) {
//            gene = "FGFR1OP-FGFR1";
//        } else if (gene.equals("PDGFRA-FIP1L1")) {
//            gene = "FIP1L1-PDGFRA";
//        } else if (gene.equals("PDGFB-COL1A1")) {
//            gene = "COL1A1-PDGFB";
//        } else if (gene.equals("BRD4-C15orf55")) {
//            gene = "BRD4-NUTM1";
//        } else if (gene.equals("MLL-MLLT3")) {
//            gene = "KMT2A-MLLT3";
//        } else if (gene.equals("BCR-ABL")) {
//            gene = "BCR-ABL1";
//        } else if (gene.equals("ERLIN2?FGFR1")) {
//            gene = "ERLIN2-FGFR1";
//        } else if (gene.equals("FGFR2?PPHLN1")) {
//            gene = "FGFR2-PPHLN1";
//        } else if (gene.equals("FGFR3 - BAIAP2L1")) {
//            gene = "FGFR3-BAIAP2L1";
//        } else if (gene.equals("PAX8-PPAR?")) {
//            gene = "PAX8-PPARA";
//        } else if (gene.contains("IGH") || gene.contains("IGK") || gene.contains("TRB") || gene.contains("Delta")
//                || gene.equals("RET-TPCN1") || gene.equals("PVT1-MYC") || gene.equals("ESR1-CCDC170") || gene.equals("BRAF-CUL1")
//                || function.contains("Loss-of-function")) {
//            gene = Strings.EMPTY;
//            typeEvent = Strings.EMPTY;
//        }
//        return ImmutableKnownFusions.builder()
//                .gene(gene)
//                .eventType(typeEvent)
//                .source(source.toString())
//                .sourceLink(source.toString())
//                .build();
//    }
}
