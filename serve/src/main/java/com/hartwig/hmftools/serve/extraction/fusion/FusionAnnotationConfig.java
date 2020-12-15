package com.hartwig.hmftools.serve.extraction.fusion;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class FusionAnnotationConfig {

    public static final Map<String, KnownFusionPair> EXONIC_FUSIONS_MAP = exonicFusionMap();

    public static final Map<String, KnownFusionPair> ODDLY_NAMED_GENES_MAP = oddlyNamedGenesMap();

    @NotNull
    private static Map<String, KnownFusionPair> oddlyNamedGenesMap() {
        Map<String, KnownFusionPair> map = Maps.newHashMap();
        KnownFusionPair fusionIGH_NKX2_1 = ImmutableKnownFusionPair.builder().geneUp("IGH").geneDown("NKX2-1").build();

        map.put("IGH-NKX2-1 Fusion", fusionIGH_NKX2_1);

        return map;
    }

    @NotNull
    private static Map<String, KnownFusionPair> exonicFusionMap() {
        Map<String, KnownFusionPair> map = Maps.newHashMap();

        KnownFusionPair fusionEGFRKDD = ImmutableKnownFusionPair.builder()
                .geneUp("EGFR")
                .minExonUp(25)
                .maxExonUp(26)
                .geneDown("EGFR")
                .minExonDown(14)
                .maxExonDown(18)
                .build();
        map.put("EGFR-KDD", fusionEGFRKDD);

        KnownFusionPair fusionEGFRvII = ImmutableKnownFusionPair.builder()
                .geneUp("EGFR")
                .minExonUp(13)
                .maxExonUp(13)
                .geneDown("EGFR")
                .minExonDown(16)
                .maxExonDown(16)
                .build();
        map.put("EGFRvII", fusionEGFRvII);

        KnownFusionPair fusionEGFRvIII = ImmutableKnownFusionPair.builder()
                .geneUp("EGFR")
                .minExonUp(1)
                .maxExonUp(1)
                .geneDown("EGFR")
                .minExonDown(8)
                .maxExonDown(8)
                .build();
        map.put("EGFRvIII", fusionEGFRvIII);
        map.put("VIII", fusionEGFRvIII);

        KnownFusionPair fusionEGFRvV = ImmutableKnownFusionPair.builder()
                .geneUp("EGFR")
                .minExonUp(24)
                .maxExonUp(24)
                .geneDown("EGFR")
                .minExonDown(29)
                .maxExonDown(29)
                .build();
        map.put("EGFRvV", fusionEGFRvV);

        KnownFusionPair fusionKITExon11 = ImmutableKnownFusionPair.builder()
                .geneUp("KIT")
                .minExonUp(11)
                .maxExonUp(11)
                .geneDown("KIT")
                .minExonDown(11)
                .maxExonDown(11)
                .build();
        map.put("EXON 11 MUTATION", fusionKITExon11);
        map.put("Exon 11 mutations", fusionKITExon11);
        map.put("Exon 11 deletions", fusionKITExon11);

        KnownFusionPair fusionMETExon14 = ImmutableKnownFusionPair.builder()
                .geneUp("MET")
                .minExonUp(13)
                .maxExonUp(13)
                .geneDown("MET")
                .minExonDown(15)
                .maxExonDown(15)
                .build();
        map.put("EXON 14 SKIPPING MUTATION", fusionMETExon14);

        return map;
    }

    private FusionAnnotationConfig() {
    }
}
