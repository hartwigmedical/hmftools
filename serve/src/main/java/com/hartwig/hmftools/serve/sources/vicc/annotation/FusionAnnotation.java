package com.hartwig.hmftools.serve.sources.vicc.annotation;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.serve.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;

import org.jetbrains.annotations.NotNull;

public final class FusionAnnotation {

    @NotNull
    public static Map<String, KnownFusionPair> createFusionExonMap() {
        Map<String, KnownFusionPair> fusionExonMap = Maps.newHashMap();

        KnownFusionPair fusionEGFR_KDD =
                ImmutableKnownFusionPair.builder().geneUp("EGFR").minExonUp(25).maxExonUp(26).geneDown("EGFR").minExonDown(14).maxExonDown(18).build();
        fusionExonMap.put("EGFR-KDD", fusionEGFR_KDD);

        KnownFusionPair fusionEGFRvII =
                ImmutableKnownFusionPair.builder().geneUp("EGFR").minExonUp(13).maxExonUp(13).geneDown("EGFR").minExonDown(16).maxExonDown(16).build();
        fusionExonMap.put("EGFRvII", fusionEGFRvII);

        KnownFusionPair fusionEGFRvV =
                ImmutableKnownFusionPair.builder().geneUp("EGFR").minExonUp(24).maxExonUp(24).geneDown("EGFR").minExonDown(29).maxExonDown(29).build();
        fusionExonMap.put("EGFRvV", fusionEGFRvV);

        KnownFusionPair fusionEGFRVIII =
                ImmutableKnownFusionPair.builder().geneUp("EGFR").minExonUp(1).maxExonUp(1).geneDown("EGFR").minExonDown(8).maxExonDown(8).build();
        fusionExonMap.put("VIII", fusionEGFRVIII);

        KnownFusionPair fusionEGFRvIII =
                ImmutableKnownFusionPair.builder().geneUp("EGFR").minExonUp(1).maxExonUp(1).geneDown("EGFR").minExonDown(8).maxExonDown(8).build();
        fusionExonMap.put("EGFRvIII", fusionEGFRvIII);

        KnownFusionPair fusionTRB_NKX2_1 =
                ImmutableKnownFusionPair.builder().geneUp("TRB").geneDown("NKX2-1").build();
        fusionExonMap.put("TRB-NKX2-1 Fusion", fusionTRB_NKX2_1);

        KnownFusionPair fusionKIT_exon11MUT =
                ImmutableKnownFusionPair.builder().geneUp("KIT").minExonUp(11).maxExonUp(11).geneDown("KIT").minExonDown(11).maxExonDown(11).build();
        fusionExonMap.put("KIT EXON 11 MUTATION", fusionKIT_exon11MUT);

        KnownFusionPair fusionKIT_exon11mut =
                ImmutableKnownFusionPair.builder().geneUp("KIT").minExonUp(11).maxExonUp(11).geneDown("KIT").minExonDown(11).maxExonDown(11).build();
        fusionExonMap.put("KIT Exon 11 mutations", fusionKIT_exon11mut);

        KnownFusionPair fusionKIT_exon11del =
                ImmutableKnownFusionPair.builder().geneUp("KIT").minExonUp(11).maxExonUp(11).geneDown("KIT").minExonDown(11).maxExonDown(11).build();
        fusionExonMap.put("KIT Exon 11 deletions", fusionKIT_exon11del);

        KnownFusionPair fusionMETexon14 =
                ImmutableKnownFusionPair.builder().geneUp("MET").minExonUp(13).maxExonUp(13).geneDown("MET").minExonDown(13).maxExonDown(13).build();
        fusionExonMap.put("MET EXON 14 SKIPPING MUTATION", fusionMETexon14);

        return fusionExonMap;
    }
}
