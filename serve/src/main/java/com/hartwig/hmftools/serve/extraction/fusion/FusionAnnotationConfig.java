package com.hartwig.hmftools.serve.extraction.fusion;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public final class FusionAnnotationConfig {

    public static final Map<String, KnownFusionPair> EXONIC_FUSIONS_MAP = exonicFusionMap();

    public static final Map<String, ExonicDelDupType> DEL_DUP_TYPE_PER_GENE = exonicDelDupTypes();

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

        KnownFusionPair fusionEGFRKinaseDomain = ImmutableKnownFusionPair.builder()
                .geneUp("EGFR")
                .minExonUp(25)
                .maxExonUp(26)
                .geneDown("EGFR")
                .minExonDown(14)
                .maxExonDown(18)
                .build();
        map.put("KINASE DOMAIN DUPLICATION (EXON 18-25)", fusionEGFRKinaseDomain);

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

        return map;
    }

    @NotNull
    private static Map<String, ExonicDelDupType> exonicDelDupTypes() {
        Map<String, ExonicDelDupType> exonicDelDupTypeMap = Maps.newHashMap();
        exonicDelDupTypeMap.put("KIT", ExonicDelDupType.PARTIAL_EXONIC_DELETION);
        exonicDelDupTypeMap.put("MET", ExonicDelDupType.FULL_EXONIC_DELETION);
        exonicDelDupTypeMap.put("EGFR", ExonicDelDupType.FULL_EXONIC_DELETION);
        exonicDelDupTypeMap.put("CBL", ExonicDelDupType.FULL_EXONIC_DELETION);
        return exonicDelDupTypeMap;
    }

    private FusionAnnotationConfig() {
    }
}
