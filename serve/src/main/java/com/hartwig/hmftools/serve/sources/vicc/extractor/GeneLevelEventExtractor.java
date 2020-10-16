package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;

import com.google.common.collect.Maps;
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
    public GeneLevelEvent extractGeneLevelEvent(@NotNull String event, @NotNull Feature feature) {
        if (event.equals("promiscuousFusion")) {
            return GeneLevelEvent.FUSION;
        } else if (event.equals("geneLevel")) {
            // TODO Determine ACTIVATION/INACTIVATION

            return GeneLevelEvent.ACTIVATION;
        } else {
            LOGGER.warn("No gene level event could be extracted from event {}", feature);
            return GeneLevelEvent.UNKONWN;
        }
    }

    @NotNull
    public Map<Feature, GeneLevelAnnotation> extractKnownGeneLevelEvents(@NotNull ViccEntry viccEntry) {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.GENE_LEVEL) {
                geneLevelEventsPerFeature.put(feature,
                        ImmutableGeneLevelAnnotation.builder().gene(feature.geneSymbol()).event(extractGeneLevelEvent("geneLevel", feature)).build());

            } else if (feature.type() == FeatureType.FUSION_PROMISCUOUS) {
                String curatedFusion = FusionCuration.curatedFusions(feature.name());
                //TODO: check if this is needed
               // if (function.equals("Likely Loss-of-function")) {
                    //            gene = Strings.EMPTY;
                    //            typeEvent = Strings.EMPTY;
                    //        }
                geneLevelEventsPerFeature.put(feature,
                        ImmutableGeneLevelAnnotation.builder().gene(curatedFusion).event(extractGeneLevelEvent("promiscuousFusion", feature)).build());
            }

        }
        return geneLevelEventsPerFeature;
    }

}
