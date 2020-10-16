package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class SignaturesExtractor {

    @NotNull
    private SignatureName extractSignatureName(@NotNull Feature feature) {
        //TODO: extent with more signatures
        if (feature.name().equals("Microsatellite Instability-High")) {
            return SignatureName.MICROSATELLITE_UNSTABLE;
        } else {
            return SignatureName.UNKNOWN;
        }
    }

    @NotNull
    public Map<Feature, SignatureName> extractSignatures(@NotNull ViccEntry viccEntry) {
        Map<Feature, SignatureName> signaturesPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.SIGNATURE) {
                SignatureName signatureName = extractSignatureName(feature);
                signaturesPerFeature.put(feature, signatureName);
            }
        }
        return signaturesPerFeature;
    }
}
