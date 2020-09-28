package com.hartwig.hmftools.serve.vicc.signatures;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.util.FeatureType;
import com.hartwig.hmftools.vicc.util.FeatureTypeExtractor;

import org.jetbrains.annotations.NotNull;

public class SignaturesExtractor {

    @NotNull
    public Map<Feature, String> extractSignatures(@NotNull ViccEntry viccEntry) {
        Map<Feature, String> signaturesPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            FeatureType featureType = FeatureTypeExtractor.extractType(feature);
            if (featureType == FeatureType.SIGNATURE) {
                signaturesPerFeature.put(feature, feature.name());
            }
        }
        return signaturesPerFeature;
    }
}
