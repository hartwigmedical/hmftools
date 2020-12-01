package com.hartwig.hmftools.serve.sources.vicc.annotation;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.serve.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;

import org.jetbrains.annotations.NotNull;

public final class FusionAnnotationConfig {

    public static final Map<String, KnownFusionPair> EXONIC_FUSIONS_MAP = exonicFusionMap();

    @NotNull
    private static Map<String, KnownFusionPair> exonicFusionMap() {
        Map<String, KnownFusionPair> fusionExonMap = Maps.newHashMap();

        KnownFusionPair fusionEGFRKDD = ImmutableKnownFusionPair.builder()
                .geneUp("EGFR")
                .minExonUp(25)
                .maxExonUp(26)
                .geneDown("EGFR")
                .minExonDown(14)
                .maxExonDown(18)
                .build();
        fusionExonMap.put("EGFR-KDD", fusionEGFRKDD);

        KnownFusionPair fusionEGFRvII = ImmutableKnownFusionPair.builder()
                .geneUp("EGFR")
                .minExonUp(13)
                .maxExonUp(13)
                .geneDown("EGFR")
                .minExonDown(16)
                .maxExonDown(16)
                .build();
        fusionExonMap.put("EGFRvII", fusionEGFRvII);

        KnownFusionPair fusionEGFRvIII = ImmutableKnownFusionPair.builder()
                .geneUp("EGFR")
                .minExonUp(1)
                .maxExonUp(1)
                .geneDown("EGFR")
                .minExonDown(8)
                .maxExonDown(8)
                .build();
        fusionExonMap.put("EGFRvIII", fusionEGFRvIII);
        fusionExonMap.put("VIII", fusionEGFRvIII);

        KnownFusionPair fusionEGFRvV = ImmutableKnownFusionPair.builder()
                .geneUp("EGFR")
                .minExonUp(24)
                .maxExonUp(24)
                .geneDown("EGFR")
                .minExonDown(29)
                .maxExonDown(29)
                .build();
        fusionExonMap.put("EGFRvV", fusionEGFRvV);

        KnownFusionPair fusionKITExon11 = ImmutableKnownFusionPair.builder()
                .geneUp("KIT")
                .minExonUp(11)
                .maxExonUp(11)
                .geneDown("KIT")
                .minExonDown(11)
                .maxExonDown(11)
                .build();
        fusionExonMap.put("EXON 11 MUTATION", fusionKITExon11);
        fusionExonMap.put("Exon 11 mutations", fusionKITExon11);
        fusionExonMap.put("Exon 11 deletions", fusionKITExon11);

        KnownFusionPair fusionMETExon14 = ImmutableKnownFusionPair.builder()
                .geneUp("MET")
                .minExonUp(13)
                .maxExonUp(13)
                .geneDown("MET")
                .minExonDown(15)
                .maxExonDown(15)
                .build();
        fusionExonMap.put("EXON 14 SKIPPING MUTATION", fusionMETExon14);

        return fusionExonMap;
    }
}
