package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.serve.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.sources.vicc.curation.FusionCuration;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class FusionExtractor {
    private static final Logger LOGGER = LogManager.getLogger(FusionExtractor.class);

    public FusionExtractor() {
    }

    public Map<Feature, KnownFusionPair> extractKnownFusions(@NotNull ViccEntry viccEntry) {
        Map<Feature, KnownFusionPair> fusionsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            String fusion = feature.name();
            String fusionGeneStart = Strings.EMPTY;
            String fusionGeneEnd = Strings.EMPTY;
            Integer exonUp = null;
            Integer exonDown = null;
            if (feature.type() == FeatureType.FUSION_PAIR) {

                if (fusion.contains("fusion")) {
                    fusion = fusion.substring(0, fusion.lastIndexOf(" "));
                }
                String curatedFusion = FusionCuration.curatedFusions(fusion);
                String[] fusionArray = curatedFusion.split("-");

                if (fusionArray.length == 2) {
                    fusionGeneStart = fusionArray[0];
                    fusionGeneEnd = fusionArray[1];
                } else if (fusionArray.length == 1) {
                    if (curatedFusion.equals("EGFRvII")) { //TODO implement correct genes and exons
                        fusionGeneStart = "EGFR";
                        exonUp = 0;
                        fusionGeneEnd = "EGFR";
                        exonDown = 0;
                    } else if (curatedFusion.equals("EGFRvV")) { //TODO implement correct genes and exons
                        fusionGeneStart = "EGFR";
                        exonUp = 0;
                        fusionGeneEnd = "EGFR";
                        exonDown = 0;
                    } else if (curatedFusion.equals("EGFRvIII") || curatedFusion.equals("VIII")) {
                        fusionGeneStart = "EGFR";
                        exonUp = 1;
                        fusionGeneEnd = "EGFR";
                        exonDown = 8;
                    } else if (curatedFusion.equals("ITD")) { //TODO implement correct genes and exons
                        fusionGeneStart = "ITD";
                        fusionGeneEnd = "ITD";
                    } else {
                        LOGGER.warn("Fusion {] is not yet curated!", curatedFusion);
                    }
                } else {
                    LOGGER.warn("Fusion name {} is to long!", curatedFusion);
                }

                fusionsPerFeature.put(feature,
                        ImmutableKnownFusionPair.builder()
                                .geneUp(fusionGeneStart)
                                .exonUp(exonUp)
                                .geneDown(fusionGeneEnd)
                                .exonDown(exonDown)
                                .build());
            } else if (feature.type() == FeatureType.FUSION_PAIR_AND_GENE_RANGE_EXON) {
                if (feature.description().equals("KIT EXON 11 MUTATION") || feature.description().equals("KIT Exon 11 mutations") || feature
                        .description()
                        .equals("KIT Exon 11 deletions")) {
                    fusionGeneStart = "KIT";
                    exonUp = 11;
                    fusionGeneEnd = "KIT";
                    exonDown = 11;
                } else if (feature.description().equals("MET EXON 14 SKIPPING MUTATION")) {
                    fusionGeneStart = "MET";
                    exonUp = 13;
                    fusionGeneEnd = "MET";
                    exonDown = 15;
                }
                fusionsPerFeature.put(feature,
                        ImmutableKnownFusionPair.builder()
                                .geneUp(fusionGeneStart)
                                .exonUp(exonUp)
                                .geneDown(fusionGeneEnd)
                                .exonDown(exonDown)
                                .build());
            }
        }
        return fusionsPerFeature;
    }

}
