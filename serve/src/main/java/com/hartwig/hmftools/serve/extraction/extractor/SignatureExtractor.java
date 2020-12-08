package com.hartwig.hmftools.serve.extraction.extractor;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SignatureExtractor {

    private static final Logger LOGGER = LogManager.getLogger(SignatureExtractor.class);

    @NotNull
    private final Set<String> microsatelliteUnstableEvents;
    @NotNull
    private final Set<String> highTumorMutationalLoadEvents;
    @NotNull
    private final Set<String> hrDeficiencyEvents;

    public SignatureExtractor(@NotNull final Set<String> microsatelliteUnstableEvents,
            @NotNull final Set<String> highTumorMutationalLoadEvents, @NotNull final Set<String> hrDeficiencyEvents) {
        this.microsatelliteUnstableEvents = microsatelliteUnstableEvents;
        this.highTumorMutationalLoadEvents = highTumorMutationalLoadEvents;
        this.hrDeficiencyEvents = hrDeficiencyEvents;
    }

    @Nullable
    public SignatureName extract(@NotNull EventType type, @NotNull String event) {
        if (type == EventType.SIGNATURE) {
            SignatureName signature = determineSignature(event);
            if (signature == null) {
                LOGGER.warn("Could not extract signature from '{}'", event);
            }
            return signature;
        }

        return null;
    }

    @Nullable
    @VisibleForTesting
    SignatureName determineSignature(@NotNull String event) {
        if (microsatelliteUnstableEvents.contains(event)) {
            return SignatureName.MICROSATELLITE_UNSTABLE;
        } else if (highTumorMutationalLoadEvents.contains(event)) {
            return SignatureName.HIGH_TUMOR_MUTATIONAL_LOAD;
        } else if (hrDeficiencyEvents.contains(event)) {
            return SignatureName.HOMOLOGOUS_RECOMBINATION_DEFICIENCY;
        }

        return null;
    }
}
