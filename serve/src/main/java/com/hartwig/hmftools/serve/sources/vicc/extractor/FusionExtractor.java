package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
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
    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    private static final Logger LOGGER = LogManager.getLogger(FusionExtractor.class);

    public FusionExtractor(@NotNull Map<String, HmfTranscriptRegion> transcriptPerGeneMap) {
        this.transcriptPerGeneMap = transcriptPerGeneMap;
    }

    public Map<Feature, KnownFusionPair> extractKnownFusions(@NotNull ViccEntry viccEntry) {
        Map<Feature, KnownFusionPair> fusionsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            String fusion = feature.name();
            String fusionGeneStart = Strings.EMPTY;
            String fusionGeneEnd = Strings.EMPTY;
            Integer minExonUp = null;
            Integer maxExonUp = null;
            Integer minExonDown = null;
            Integer maxExonDown = null;

            if (feature.type() == FeatureType.FUSION_PAIR) {
                String[] fusionArray = fusion.split("-");

                if (fusionArray.length == 2) {
                    if (fusion.equals("EGFR-KDD")) {
                        fusionGeneStart = feature.geneSymbol();
                        minExonUp = 25;
                        maxExonUp = 26;
                        fusionGeneEnd = feature.geneSymbol();
                        minExonDown = 14;
                        maxExonDown = 18;
                    } else {
                        fusionGeneStart = fusionArray[0];
                        fusionGeneEnd = fusionArray[1].split(" ")[0];
                    }
                } else if (fusionArray.length == 1) {
                    if (fusion.equals("EGFRvII")) {
                        fusionGeneStart = feature.geneSymbol();
                        minExonUp = 13;
                        maxExonUp =13;
                        fusionGeneEnd = feature.geneSymbol();
                        minExonDown = 16;
                        maxExonDown = 16;
                    } else if (fusion.equals("EGFRvV")) {
                        fusionGeneStart = feature.geneSymbol();
                        minExonUp = 24;
                        maxExonUp =24;
                        fusionGeneEnd = feature.geneSymbol();
                        minExonDown = 29;
                        maxExonDown =29;
                    } else if (fusion.equals("EGFRvIII") || fusion.equals("VIII")) {
                        fusionGeneStart = feature.geneSymbol();
                        minExonUp = 1;
                        maxExonUp =1;
                        fusionGeneEnd = feature.geneSymbol();
                        minExonDown = 8;
                        maxExonDown =8;
                    } else {
                        LOGGER.warn("Fusion '{}' can not be interpreted!", fusion);
                    }
                } else {
                    if (fusion.equals("TRB-NKX2-1 Fusion")) {
                        fusionGeneStart = "TRB";
                        fusionGeneEnd = "NKX2-1";
                    } else {
                        LOGGER.warn("Too many parts in fusion name: {}!", fusion);
                    }
                }

                HmfTranscriptRegion canonicalTranscriptStart = transcriptPerGeneMap.get(fusionGeneStart);
                HmfTranscriptRegion canonicalTranscriptEnd = transcriptPerGeneMap.get(fusionGeneEnd);

                if (canonicalTranscriptStart == null || canonicalTranscriptEnd == null) {
                    LOGGER.warn(
                            "Could not find fusion gene start {} or fusion gene end {} in HMF gene panel. Skipping fusion pair extraction!",
                            fusionGeneStart,
                            fusionGeneEnd);
                } else {
                    fusionsPerFeature.put(feature,
                            ImmutableKnownFusionPair.builder()
                                    .geneUp(fusionGeneStart)
                                    .minExonUp(minExonUp)
                                    .maxExonUp(maxExonUp)
                                    .geneDown(fusionGeneEnd)
                                    .minExonDown(minExonDown)
                                    .maxExonDown(maxExonDown)
                                    .build());
                }

            } else if (feature.type() == FeatureType.FUSION_PAIR_AND_GENE_RANGE_EXON) {
                if (feature.description().equals("KIT EXON 11 MUTATION") || feature.description().equals("KIT Exon 11 mutations") || feature
                        .description()
                        .equals("KIT Exon 11 deletions")) {
                    fusionGeneStart = feature.geneSymbol();
                    minExonUp = extractExonNumber(feature.name());
                    maxExonUp =extractExonNumber(feature.name());
                    fusionGeneEnd = feature.geneSymbol();
                                     minExonDown = extractExonNumber(feature.name());
                    maxExonDown =extractExonNumber(feature.name());
                } else if (feature.description().equals("MET EXON 14 SKIPPING MUTATION")) {
                    fusionGeneStart = feature.geneSymbol();
                    minExonUp = extractExonNumber(feature.name()) - 1;
                    maxExonUp =extractExonNumber(feature.name()) - 1;
                    fusionGeneEnd = feature.geneSymbol();
                    minExonDown = extractExonNumber(feature.name()) + 1;
                    maxExonDown =extractExonNumber(feature.name()) + 1;
                }

                HmfTranscriptRegion canonicalTranscriptStart = transcriptPerGeneMap.get(fusionGeneStart);
                HmfTranscriptRegion canonicalTranscriptEnd = transcriptPerGeneMap.get(fusionGeneEnd);

                if (canonicalTranscriptStart == null || canonicalTranscriptEnd == null) {
                    LOGGER.warn(
                            "Could not find fusion gene start {} or fusion gene end {} in HMF gene panel. Skipping fusion par and gene "
                                    + "range exon extraction for internal fusion!",
                            fusionGeneStart,
                            fusionGeneEnd);
                } else {
                    fusionsPerFeature.put(feature,
                            ImmutableKnownFusionPair.builder()
                                    .geneUp(fusionGeneStart)
                                    .minExonUp(minExonUp)
                                    .maxExonUp(maxExonUp)
                                    .geneDown(fusionGeneEnd)
                                    .minExonDown(minExonDown)
                                    .maxExonDown(maxExonDown)
                                    .build());
                }
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
