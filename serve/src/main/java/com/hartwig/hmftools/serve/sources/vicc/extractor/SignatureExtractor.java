package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SignatureExtractor {

    private static final Logger LOGGER = LogManager.getLogger(SignatureExtractor.class);

    public SignatureExtractor() {
    }

    @NotNull
    public Map<Feature, SignatureName> extract(@NotNull ViccEntry viccEntry) {
        Map<Feature, SignatureName> signaturesPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            if (feature.type() == MutationType.SIGNATURE) {
                SignatureName signatureName = extractSignatureName(feature.name());
                if (signatureName != null) {
                    signaturesPerFeature.put(feature, signatureName);
                } else {
                    LOGGER.warn("Could not extract signature from '{}'", feature.name());
                }
            }
        }
        return signaturesPerFeature;
    }

    @Nullable
    @VisibleForTesting
    static SignatureName extractSignatureName(@NotNull String featureName) {
        if (featureName.equals("Microsatellite Instability-High")) {
            return SignatureName.MICROSATELLITE_UNSTABLE;
        }

        return null;
    }
}
