package com.hartwig.hmftools.serve.vicc.extractor;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.vicc.annotation.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.vicc.annotation.ImmutableGeneLevelAnnotation;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class GeneLevelEventExtractor {

    public GeneLevelEventExtractor() {
    }

    @NotNull
    public Map<Feature, GeneLevelAnnotation> extractKnownGeneLevelEvents(@NotNull ViccEntry viccEntry) {
        Map<Feature, GeneLevelAnnotation> geneLevelEventsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.GENE_LEVEL) {
                // TODO Extract GeneEvent from Feature
                geneLevelEventsPerFeature.put(feature,
                        ImmutableGeneLevelAnnotation.builder().gene(feature.geneSymbol()).event(GeneLevelEvent.ACTIVATION).build());

            }

            //             else if (isValidSingleCodonRange(feature.proteinAnnotation())) {
            //                LOGGER.info(feature.proteinAnnotation());
            //                geneLevelEventsPerFeature.put(feature, feature.geneSymbol());
            //            }

            //            if (GENE_ACTIVATION.contains(feature.name()) || GENE_ACTIVATION.contains(feature.biomarkerType()) || GENE_ACTIVATION.contains(
            //                    feature.proteinAnnotation())) {
            //                geneLevelEventsPerFeature.put(feature, "gain of " + feature.geneSymbol());
            //            } else if (GENE_INACTIVATION.contains(feature.name()) || GENE_INACTIVATION.contains(feature.biomarkerType())
            //                    || GENE_INACTIVATION.contains(feature.provenanceRule()) || GENE_INACTIVATION.contains(feature.proteinAnnotation())) {
            //                geneLevelEventsPerFeature.put(feature, "loss of " + feature.geneSymbol());
            //            }
        }
        return geneLevelEventsPerFeature;
    }

}
