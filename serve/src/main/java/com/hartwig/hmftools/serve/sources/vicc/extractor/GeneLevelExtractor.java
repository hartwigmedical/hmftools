package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneLevelAnnotation;
import com.hartwig.hmftools.vicc.annotation.GeneLevelClassifier;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GeneLevelExtractor {

    private static final Logger LOGGER = LogManager.getLogger(GeneLevelExtractor.class);

    private static final String GENE_ONLY = "gene_only";

    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;
    @NotNull
    private final List<DriverGene> driverGenes;

    public GeneLevelExtractor(@NotNull final Map<String, HmfTranscriptRegion> transcriptPerGeneMap,
            @NotNull final List<DriverGene> driverGenes) {
        this.transcriptPerGeneMap = transcriptPerGeneMap;
        this.driverGenes = driverGenes;
    }

    @NotNull
    public Map<Feature, GeneLevelAnnotation> extractGeneLevelEvents(@NotNull ViccEntry viccEntry) {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(feature.geneSymbol());
            if (feature.type() == EventType.GENE_LEVEL) {
                if (canonicalTranscript == null) {
                    LOGGER.warn("Could not find gene '{}' in HMF gene panel. Skipping gene level extraction!", feature.geneSymbol());
                } else {
                    geneLevelEventsPerFeature.put(feature,
                            ImmutableGeneLevelAnnotation.builder()
                                    .gene(feature.geneSymbol())
                                    .event(extractGeneLevelEvent(feature, driverGenes))
                                    .build());
                }
            } else if (feature.type() == EventType.PROMISCUOUS_FUSION) {
                if (canonicalTranscript == null) {
                    LOGGER.warn("Could not find gene '{}' in HMF gene panel. Skipping gene level extraction!", feature.geneSymbol());
                } else {
                    geneLevelEventsPerFeature.put(feature,
                            ImmutableGeneLevelAnnotation.builder().gene(feature.geneSymbol()).event(GeneLevelEvent.FUSION).build());
                }
            }

        }

        return geneLevelEventsPerFeature;
    }

    @NotNull
    @VisibleForTesting
    static GeneLevelEvent extractGeneLevelEvent(@NotNull Feature feature, @NotNull List<DriverGene> driverGenes) {
        String eventDescription = feature.description().split(" ", 2)[1].trim();
        if (GeneLevelClassifier.GENERIC_GENE_LEVEL_KEYWORDS.contains(eventDescription) || feature.provenanceRule() != null) {
            for (DriverGene driverGene : driverGenes) {
                if (driverGene.gene().equals(feature.geneSymbol())) {
                    if (driverGene.likelihoodType() == DriverCategory.ONCO) {

                        if (feature.provenanceRule() != null) {
                            if (feature.provenanceRule().equals(GENE_ONLY)) {
                                return GeneLevelEvent.ACTIVATION;
                            } else {
                                return GeneLevelEvent.ACTIVATION;
                            }
                        } else if (GeneLevelClassifier.GENERIC_GENE_LEVEL_KEYWORDS.contains(eventDescription)) {
                            return GeneLevelEvent.ACTIVATION;
                        }
                    } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                        if (feature.provenanceRule() != null) {
                            if (feature.provenanceRule().equals(GENE_ONLY)) {
                                return GeneLevelEvent.INACTIVATION;
                            } else {
                                return GeneLevelEvent.INACTIVATION;
                            }
                        } else if (GeneLevelClassifier.GENERIC_GENE_LEVEL_KEYWORDS.contains(eventDescription)) {
                            return GeneLevelEvent.INACTIVATION;
                        }
                    }
                }
            }
            LOGGER.warn("Gene {} is not present in driver catalog in gene level extractor", feature.geneSymbol());
        } else if (GeneLevelClassifier.INACTIVATING_GENE_LEVEL_KEYWORDS.contains(eventDescription)) {
            return GeneLevelEvent.INACTIVATION;
        } else if (GeneLevelClassifier.ACTIVATING_GENE_LEVEL_KEYWORDS.contains(eventDescription)) {
            return GeneLevelEvent.ACTIVATION;
        } else {
            // LOGGER.warn("Unknown event {}", feature);
            return GeneLevelEvent.UNKNOWN;
        }
        //  LOGGER.warn("Unknown event {}", feature);
        return GeneLevelEvent.UNKNOWN;
    }
}
