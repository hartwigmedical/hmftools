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
    public static GeneLevelEvent extractGeneLevelEvent(@NotNull String event, @NotNull Feature feature,
            @NotNull List<DriverGene> driverGenes) {
        if (event.equals("promiscuousFusion")) {
            return GeneLevelEvent.FUSION;
        } else if (event.equals("geneLevel")) {
            GeneLevelEvent geneLevelEvent = GeneLevelEvent.UNKONWN;
            for (DriverGene driverGene : driverGenes) {
                if (driverGene.gene().equals(feature.geneSymbol())) {
                    if (driverGene.likelihoodType() == DriverCategory.TSG) {
                        geneLevelEvent = GeneLevelEvent.ACTIVATION;
                        // TODO Determine ACTIVATION/INACTIVATION for TSG
                    } else if (driverGene.likelihoodType() == DriverCategory.ONCO) {
                        // TODO Determine ACTIVATION/INACTIVATION for ONCO
                        geneLevelEvent = GeneLevelEvent.ACTIVATION;
                    }
                } else {
                    //TODO why is this true for values which are present in driver catalog
                   // LOGGER.warn("Gene {} is not present in driver catalog and could not determine driver category", feature.geneSymbol());
                }
            }
            return geneLevelEvent;
        } else {
            LOGGER.warn("No gene level event could be extracted from event {}", feature);
            return GeneLevelEvent.UNKONWN;
        }
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
                                .event(extractGeneLevelEvent("geneLevel", feature, driverGenes))
                                .build());

            } else if (feature.type() == FeatureType.FUSION_PROMISCUOUS) {

                String curatedPromiscuousFusion = FusionCuration.curatedFusions(feature.geneSymbol());
                //TODO: check if this is needed
                // if (function.equals("Likely Loss-of-function")) {
                //            gene = Strings.EMPTY;
                //            typeEvent = Strings.EMPTY;
                //        }
                LOGGER.info(feature);
                geneLevelEventsPerFeature.put(feature,
                        ImmutableGeneLevelAnnotation.builder()
                                .gene(curatedPromiscuousFusion)
                                .event(extractGeneLevelEvent("promiscuousFusion", feature, driverGenes))
                                .build());
            }

        }
        return geneLevelEventsPerFeature;
    }

}
