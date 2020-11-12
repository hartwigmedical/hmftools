package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.serve.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.fusion.KnownFusionPair;
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
                String[] fusionArray = fusion.split("-");

                if (fusionArray.length == 2) {
                    fusionGeneStart = fusionArray[0];
                    fusionGeneEnd = fusionArray[1].split(" ")[0];
                } else if (fusionArray.length == 1) {
                    if (fusion.equals("EGFRvII")) { //TODO implement correct genes and exons
                        fusionGeneStart = feature.geneSymbol();
                        exonUp = 0;
                        fusionGeneEnd = feature.geneSymbol();
                        exonDown = 0;
                    } else if (fusion.equals("EGFRvV")) { //TODO implement correct genes and exons
                        fusionGeneStart = feature.geneSymbol();
                        exonUp = 0;
                        fusionGeneEnd = feature.geneSymbol();
                        exonDown = 0;
                    } else if (fusion.equals("EGFRvIII") || fusion.equals("VIII")) {
                        fusionGeneStart = feature.geneSymbol();
                        exonUp = 1;
                        fusionGeneEnd = feature.geneSymbol();
                        exonDown = 8;
                    } else {
                        LOGGER.warn("Fusion '{}' can not be interpreted!", fusion);
                    }
                } else {
                    LOGGER.warn("Too many parts in fusion name: {}!", fusion);
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
                    fusionGeneStart = feature.geneSymbol();
                    exonUp = extractExonNumber(feature.name());
                    fusionGeneEnd = feature.geneSymbol();
                    exonDown = extractExonNumber(feature.name());
                } else if (feature.description().equals("MET EXON 14 SKIPPING MUTATION")) {
                    fusionGeneStart = feature.geneSymbol();
                    exonUp = extractExonNumber(feature.name()) -1;
                    fusionGeneEnd = feature.geneSymbol();
                    exonDown = extractExonNumber(feature.name()) + 1;
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

    @NotNull
    public static Integer extractExonNumber(@NotNull String featureName) {
        String exonNumberAsString = featureName.split(" ")[1];
        return Integer.valueOf(exonNumberAsString);
    }
}
