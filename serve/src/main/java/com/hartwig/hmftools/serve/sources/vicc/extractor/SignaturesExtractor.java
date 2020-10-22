package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class SignaturesExtractor {

    @NotNull
    @VisibleForTesting
    public static SignatureName extractSignatureName(@NotNull String featureName) {
        if (featureName.equals("Microsatellite Instability-High")) {
            return SignatureName.MICROSATELLITE_UNSTABLE;
        } else if (featureName.equals("High mutational load")) { // TODO: check for right string
            return SignatureName.HIGH_TUMOR_MUTATIONAL_LOAD;
        } else if (featureName.equals("HRD")) {
            return SignatureName.HOMOLOGOUS_RECOMBINATION_DEFICIENY; // TODO: check for right string
        } else {
            return SignatureName.UNKNOWN;
        }
    }

    @NotNull
    public Map<Feature, SignatureName> extractSignatures(@NotNull ViccEntry viccEntry) {
        Map<Feature, SignatureName> signaturesPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.SIGNATURE) {
                SignatureName signatureName = extractSignatureName(feature.name());
                signaturesPerFeature.put(feature, signatureName);
            }
        }
        return signaturesPerFeature;
    }
}
