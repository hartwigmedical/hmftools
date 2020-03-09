package com.hartwig.hmftools.knowledgebasegenerator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum GenomicEvents {
    AMPLIFICATION,
    DELETION,
    VARIANT,
    RANGE,
    FUSION,
    SIGNATURE,
    UNKNOWN;

    private static final Logger LOGGER = LogManager.getLogger(GenomicEvents.class);

    @NotNull
    public static GenomicEvents genomicEvents(@NotNull String typeEvent) {
        if (typeEvent.equals("Amplification")) {
            return AMPLIFICATION;
        } else if (typeEvent.equals("Deletion")) {
            return DELETION;
        } else if (typeEvent.equals("Variants")) {
            return VARIANT;
        } else if (typeEvent.equals("Range")) {
            return RANGE;
        } else if (typeEvent.equals("Fusions")) {
            return FUSION;
        } else if (typeEvent.equals("Signatures")) {
            return SIGNATURE;
        } else {
            LOGGER.warn("Unknown genomic event!");
            return UNKNOWN;
        }
    }
}
