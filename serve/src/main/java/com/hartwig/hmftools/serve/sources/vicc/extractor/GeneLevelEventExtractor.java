package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneLevelAnnotation;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.annotation.GeneRangeClassifier;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GeneLevelEventExtractor {

    private static final Logger LOGGER = LogManager.getLogger(GeneLevelEventExtractor.class);

    private static final String GENE_ONLY = "gene_only";
    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    public GeneLevelEventExtractor(@NotNull Map<String, HmfTranscriptRegion> transcriptPerGeneMap) {
        this.transcriptPerGeneMap = transcriptPerGeneMap;
    }

    @NotNull
    public static GeneLevelEvent extractGeneLevelEvent(@NotNull Feature feature, @NotNull List<DriverGene> driverGenes) {
        String eventDescription = feature.description().split(" ", 2)[1].trim();
        if (GeneRangeClassifier.GENERIC_GENE_LEVEL_KEYWORDS.contains(eventDescription) || feature.provenanceRule() != null) {
            for (DriverGene driverGene : driverGenes) {
                if (driverGene.gene().equals(feature.geneSymbol())) {
                    if (driverGene.likelihoodType() == DriverCategory.ONCO) {

                        if (feature.provenanceRule() != null) {
                            if (feature.provenanceRule().equals(GENE_ONLY)) {
                                return GeneLevelEvent.ACTIVATION;
                            } else {
                                return GeneLevelEvent.ACTIVATION;
                            }
                        } else if (GeneRangeClassifier.GENERIC_GENE_LEVEL_KEYWORDS.contains(eventDescription)) {
                            return GeneLevelEvent.ACTIVATION;
                        }
                    } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                        if (feature.provenanceRule() != null) {
                            if (feature.provenanceRule().equals(GENE_ONLY)) {
                                return GeneLevelEvent.INACTIVATION;
                            } else {
                                return GeneLevelEvent.INACTIVATION;
                            }
                        } else if (GeneRangeClassifier.GENERIC_GENE_LEVEL_KEYWORDS.contains(eventDescription)) {
                            return GeneLevelEvent.INACTIVATION;
                        }
                    }
                }
            }
            LOGGER.warn("Gene {} is not present in driver catalog in gene level extractor", feature.geneSymbol());
        } else if (GeneRangeClassifier.INACTIVATING_GENE_LEVEL_KEYWORDS.contains(eventDescription)) {
            return GeneLevelEvent.INACTIVATION;
        } else if (GeneRangeClassifier.ACTIVATING_GENE_LEVEL_KEYWORDS.contains(eventDescription)) {
            return GeneLevelEvent.ACTIVATION;
        } else {
           // LOGGER.warn("Unknown event {}", feature);
            return GeneLevelEvent.UNKNOWN;
        }
      //  LOGGER.warn("Unknown event {}", feature);
        return GeneLevelEvent.UNKNOWN;
    }

    @NotNull
    public Map<Feature, GeneLevelAnnotation> extractKnownGeneLevelEvents(@NotNull ViccEntry viccEntry,
            @NotNull List<DriverGene> driverGenes) {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(feature.geneSymbol());

            if (feature.type() == FeatureType.GENE_LEVEL) {

                if (canonicalTranscript == null) {
                    LOGGER.warn("Could not find gene {} in HMF gene panel. Skipping gene level extraction!", feature.geneSymbol());
                } else {
                    geneLevelEventsPerFeature.put(feature,
                            ImmutableGeneLevelAnnotation.builder()
                                    .gene(feature.geneSymbol())
                                    .event(extractGeneLevelEvent(feature, driverGenes))
                                    .build());
                }
            } else if (feature.type() == FeatureType.PROMISCUOUS_FUSION) {

                if (canonicalTranscript == null) {
                    LOGGER.warn("Could not find gene {} in HMF gene panel. Skipping promiscuous fusion extraction!", feature.geneSymbol());
                } else {
                    geneLevelEventsPerFeature.put(feature,
                            ImmutableGeneLevelAnnotation.builder().gene(feature.geneSymbol()).event(GeneLevelEvent.FUSION).build());
                }

            }

        }
        return geneLevelEventsPerFeature;
    }

}
