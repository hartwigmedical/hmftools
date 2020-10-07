package com.hartwig.hmftools.serve.vicc.signatures;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class SignaturesExtractor {

    @NotNull
    public Map<Feature, SignatureName> extractSignatures(@NotNull ViccEntry viccEntry) {
        Map<Feature, SignatureName> signaturesPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.SIGNATURE) {
                // TODO Map to SignatureName
//                signaturesPerFeature.put(feature, feature.name());
            }
        }
        return signaturesPerFeature;
    }
}
