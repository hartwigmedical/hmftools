package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.sources.vicc.annotation.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.annotation.ImmutableGeneLevelAnnotation;
import com.hartwig.hmftools.serve.sources.vicc.curation.FusionCuration;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.annotation.GeneRangeClassifier;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GeneLevelEventExtractor {

    private static final Logger LOGGER = LogManager.getLogger(GeneLevelEventExtractor.class);

    public GeneLevelEventExtractor() {
    }

    @NotNull
    public static GeneLevelEvent extractGeneLevelEvent(@NotNull Feature feature, @NotNull List<DriverGene> driverGenes) {
        String eventDescription = feature.description().split(" ", 2)[1].trim();
        if (GeneRangeClassifier.DETAILED_GENE_LEVEL_INFO_WITHOUT_TSG_ONCO.contains(eventDescription) || feature.provenanceRule() != null) {
            for (DriverGene driverGene : driverGenes) {
                if (driverGene.gene().equals(feature.geneSymbol())) {
                    if (driverGene.likelihoodType() == DriverCategory.ONCO) {

                        if (feature.provenanceRule() != null) {
                            if (feature.provenanceRule().equals(GeneRangeClassifier.GENE_LEVEL)) {
                                return GeneLevelEvent.ACTIVATION;
                            } else {
                                return GeneLevelEvent.ACTIVATION;
                            }
                        } else if (GeneRangeClassifier.DETAILED_GENE_LEVEL_INFO_WITHOUT_TSG_ONCO.contains(eventDescription)) {
                            return GeneLevelEvent.ACTIVATION;
                        }
                    } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                        if (feature.provenanceRule() != null) {
                            if (feature.provenanceRule().equals(GeneRangeClassifier.GENE_LEVEL)) {
                                return GeneLevelEvent.INACTIVATION;
                            } else {
                                return GeneLevelEvent.INACTIVATION;
                            }
                        } else if (GeneRangeClassifier.DETAILED_GENE_LEVEL_INFO_WITHOUT_TSG_ONCO.contains(eventDescription)) {
                            return GeneLevelEvent.INACTIVATION;
                        }
                    }
                }
            }
            LOGGER.warn("Gene {} is not present in driver catalog", feature.geneSymbol());
        } else if (GeneRangeClassifier.DETAILED_GENE_LEVEL_INFO_WITH_TSG.contains(eventDescription)) {
            return GeneLevelEvent.INACTIVATION;
        } else if (GeneRangeClassifier.DETAILED_GENE_LEVEL_INFO_WITH_ONCO.contains(eventDescription)) {
            return GeneLevelEvent.ACTIVATION;
        } else {
            LOGGER.warn("Unknown event {}", feature);
            return GeneLevelEvent.UNKNOWN;
        }
        LOGGER.warn("Unknown event {}", feature);
        return GeneLevelEvent.UNKNOWN;
    }

    @NotNull
    public Map<Feature, GeneLevelAnnotation> extractKnownGeneLevelEvents(@NotNull ViccEntry viccEntry,
            @NotNull List<DriverGene> driverGenes) {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.GENE_LEVEL) {
                geneLevelEventsPerFeature.put(feature,
                        ImmutableGeneLevelAnnotation.builder()
                                .gene(feature.geneSymbol())
                                .event(extractGeneLevelEvent(feature, driverGenes))
                                .build());

            } else if (feature.type() == FeatureType.FUSION_PROMISCUOUS) {

                String curatedPromiscuousFusion = FusionCuration.curatedFusions(feature.geneSymbol(), feature);
                //TODO: check if this is needed
                // if (function.equals("Likely Loss-of-function")) {
                //            gene = Strings.EMPTY;
                //            typeEvent = Strings.EMPTY;
                //        }
                geneLevelEventsPerFeature.put(feature,
                        ImmutableGeneLevelAnnotation.builder().gene(curatedPromiscuousFusion).event(GeneLevelEvent.FUSION).build());
            }

        }
        return geneLevelEventsPerFeature;
    }

}
