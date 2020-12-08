package com.hartwig.hmftools.serve.sources.vicc.extractor;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SignatureExtractor {

    private static final Logger LOGGER = LogManager.getLogger(SignatureExtractor.class);

    public SignatureExtractor() {
    }

    @Nullable
    public SignatureName extract(@NotNull EventType type, @NotNull String event) {
        if (type == EventType.SIGNATURE) {
            SignatureName signatureName = extractSignatureName(event);
            if (signatureName != null) {
                return signatureName;
            } else {
                LOGGER.warn("Could not extract signature from '{}'", event);
            }
        }

        return null;
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
