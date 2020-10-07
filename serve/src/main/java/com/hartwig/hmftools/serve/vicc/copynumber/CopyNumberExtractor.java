package com.hartwig.hmftools.serve.vicc.copynumber;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class CopyNumberExtractor {

    public CopyNumberExtractor() {
    }

    @NotNull
    public Map<Feature, CopyNumberAnnotation> extractKnownAmplificationsDeletions(@NotNull ViccEntry viccEntry) {
        Map<Feature, CopyNumberAnnotation> ampsDelsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.AMPLIFICATION) {
                ampsDelsPerFeature.put(feature,
                        ImmutableCopyNumberAnnotation.builder().gene(feature.geneSymbol()).type(CopyNumberType.AMPLIFICATION).build());
            } else if (feature.type() == FeatureType.DELETION) {
                ampsDelsPerFeature.put(feature,
                        ImmutableCopyNumberAnnotation.builder().gene(feature.geneSymbol()).type(CopyNumberType.DELETION).build());
            }
        } return ampsDelsPerFeature;
    }
}
