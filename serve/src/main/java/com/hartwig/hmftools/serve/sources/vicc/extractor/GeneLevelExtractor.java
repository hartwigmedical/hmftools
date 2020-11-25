package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.common.serve.classification.matchers.GeneLevelMatcher;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.check.CheckGenes;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
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
                    CheckGenes.checkGensInPanel(feature.geneSymbol(), feature.name());
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
        String event = Strings.EMPTY;
        if (feature.name().split(" ").length > 1) {
            event = feature.name().split(" ", 2)[1].trim();
        } else {
            event = feature.name();
        }

        if (GeneLevelMatcher.INACTIVATING_GENE_LEVEL_KEYWORDS.contains(event)) {
            return GeneLevelEvent.INACTIVATION;
        } else if (GeneLevelMatcher.ACTIVATING_GENE_LEVEL_KEYWORDS.contains(event)) {
            return GeneLevelEvent.ACTIVATION;
        } else if (GeneLevelMatcher.GENERIC_GENE_LEVEL_KEYWORDS.contains(event)) {
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
        } else if (feature.provenanceRule() != null){
            if (GENE_ONLY.contains(feature.provenanceRule())) {
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
            }
        }else {
            // LOGGER.warn("Unknown event {}", feature);
            return GeneLevelEvent.UNKNOWN;
        }
        //  LOGGER.warn("Unknown event {}", feature);
        return GeneLevelEvent.UNKNOWN;
    }
}
