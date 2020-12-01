package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.check.CheckGenes;
import com.hartwig.hmftools.vicc.annotation.ViccClassificationConfig;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GeneLevelExtractor {

    private static final Logger LOGGER = LogManager.getLogger(GeneLevelExtractor.class);

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
        boolean usingGenes;

        for (Feature feature : viccEntry.features()) {
            HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(feature.geneSymbol());
            if (feature.type() == MutationType.GENE_LEVEL) {
                if (canonicalTranscript == null) {
                    CheckGenes.checkGensInPanel(feature.geneSymbol(), feature.name());

                } else {
                    geneLevelEventsPerFeature.put(feature,
                            ImmutableGeneLevelAnnotation.builder()
                                    .gene(feature.geneSymbol())
                                    .event(extractGeneLevelEvent(feature, driverGenes))
                                    .build());
                }
            } else if (feature.type() == MutationType.PROMISCUOUS_FUSION) {
                if (canonicalTranscript == null) {
                    usingGenes = CheckGenes.checkGensInPanelForCuration(feature.geneSymbol(), feature.name());
                    if (usingGenes) {
                        geneLevelEventsPerFeature.put(feature,
                                ImmutableGeneLevelAnnotation.builder().gene(feature.geneSymbol()).event(GeneLevelEvent.FUSION).build());
                    }
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
    public static GeneLevelEvent extractGeneLevelEvent(@NotNull Feature feature, @NotNull List<DriverGene> driverGenes) {
        String event;
        String geneSymbol = feature.geneSymbol();
        String geneSymbolEvent = feature.name().split(" ")[0];

        if (geneSymbolEvent.equals(geneSymbol)) {
            event = feature.name().split(" ", 2)[1].trim();

        } else {
            event = feature.name();
        }

        if (ViccClassificationConfig.INACTIVATING_GENE_LEVEL_KEY_PHRASES.contains(event)) {
            return GeneLevelEvent.INACTIVATION;
        } else if (ViccClassificationConfig.ACTIVATING_GENE_LEVEL_KEY_PHRASES.contains(event)) {
            return GeneLevelEvent.ACTIVATION;
        } else if (ViccClassificationConfig.GENERIC_GENE_LEVEL_KEY_PHRASES.contains(event)) {
            return extractGeneLevelEventGene(feature, driverGenes);
        } else if (feature.geneSymbol().equals(feature.name().replaceAll("\\s+", ""))) {
            return extractGeneLevelEventGene(feature, driverGenes);
        } else {
            LOGGER.warn("Unknown event {}", feature);
            return GeneLevelEvent.UNKNOWN;
        }
    }

    @VisibleForTesting
    @NotNull
    public static GeneLevelEvent extractGeneLevelEventGene(@NotNull Feature feature, @NotNull List<DriverGene> driverGenes) {
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(feature.geneSymbol())) {
                if (driverGene.likelihoodType() == DriverCategory.ONCO) {
                    return GeneLevelEvent.ACTIVATION;
                } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                    return GeneLevelEvent.INACTIVATION;
                }
            }
        }
        CheckGenes.checkGensInPanel(feature.geneSymbol(), feature.name());
        return GeneLevelEvent.UNKNOWN;
    }
}
