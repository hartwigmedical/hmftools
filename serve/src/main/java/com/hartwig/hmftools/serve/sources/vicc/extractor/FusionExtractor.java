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

    @NotNull
    public Map<Feature, KnownFusionPair> extractFusionPairs(@NotNull ViccEntry viccEntry,
            @NotNull Map<String, KnownFusionPair> fusionAnnotations) {
        Map<Feature, KnownFusionPair> fusionsPerFeature = Maps.newHashMap();
        KnownFusionPair annotatedFusion = ImmutableKnownFusionPair.builder().geneUp(Strings.EMPTY).geneDown(Strings.EMPTY).build();

        for (Feature feature : viccEntry.features()) {
            String fusion = feature.name();

            if (feature.type() == FeatureType.FUSION_PAIR) {
                String[] fusionArray = fusion.split("-");

                if (fusionArray.length == 2) {
                    if (fusionAnnotations.containsKey(fusion)) {
                        annotatedFusion = ImmutableKnownFusionPair.builder().from(fusionAnnotations.get(fusion)).build();

                    } else {
                        annotatedFusion = ImmutableKnownFusionPair.builder().geneUp(fusionArray[0]).geneDown(fusionArray[1].split(" ")[0]).build();
                    }
                } else if (fusionArray.length == 1) {
                    if (fusionAnnotations.containsKey(fusion)) {
                        annotatedFusion = ImmutableKnownFusionPair.builder().from(fusionAnnotations.get(fusion)).build();
                    } else {
                        LOGGER.warn("Fusion '{}' can not be interpreted!", fusion);
                    }
                } else {
                    if (fusionAnnotations.containsKey(fusion)) {
                        annotatedFusion = ImmutableKnownFusionPair.builder().from(fusionAnnotations.get(fusion)).build();
                    } else {
                        LOGGER.warn("Too many parts in fusion name: {}!", fusion);
                    }
                }

                HmfTranscriptRegion canonicalTranscriptStart = transcriptPerGeneMap.get(annotatedFusion.geneUp());
                HmfTranscriptRegion canonicalTranscriptEnd = transcriptPerGeneMap.get(annotatedFusion.geneDown());

                if (canonicalTranscriptStart == null || canonicalTranscriptEnd == null) {
                    LOGGER.warn(
                            "Could not find fusion gene start {} or fusion gene end {} in HMF gene panel. Skipping fusion pair extraction!",
                            annotatedFusion.geneUp(),
                            annotatedFusion.geneDown());
                } else {
                    fusionsPerFeature.put(feature, annotatedFusion);
                }

            } else if (feature.type() == FeatureType.FUSION_PAIR_AND_GENE_RANGE_EXON) {
                if (fusionAnnotations.containsKey(feature.description())) {
                    annotatedFusion = ImmutableKnownFusionPair.builder().from(fusionAnnotations.get(feature.description())).build();
                }

                HmfTranscriptRegion canonicalTranscriptStart = transcriptPerGeneMap.get(annotatedFusion.geneUp());
                HmfTranscriptRegion canonicalTranscriptEnd = transcriptPerGeneMap.get(annotatedFusion.geneDown());

                if (canonicalTranscriptStart == null || canonicalTranscriptEnd == null) {
                    LOGGER.warn("Could not find fusion gene start {} or fusion gene end {} in HMF gene panel. Skipping fusion par and gene "
                            + "range exon extraction for internal fusion!", annotatedFusion.geneUp(), annotatedFusion.geneDown());
                } else {
                    fusionsPerFeature.put(feature, annotatedFusion);
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
