package com.hartwig.hmftools.knowledgebasegenerator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum GenomicEvents {
    AMPLIFICATION,
    DELETION,
    VARIANT,
    RANGE,
    FUSION_PAIR,
    FUSION_PROMISCUOUS,
    SIGNATURE,
    SIGNATURE_MSI,
    SIGNATURE_HRD,
    SIGNATURE_MTL,
    SIGNATURE_MTB,
    UNKNOWN;

    private static final Logger LOGGER = LogManager.getLogger(GenomicEvents.class);

    @NotNull
    public static GenomicEvents genomicEvents(@NotNull String typeEvent) {
        switch (typeEvent) {
            case "Amplification":
                return AMPLIFICATION;
            case "Deletion":
                return DELETION;
            case "Variants":
                return VARIANT;
            case "Range":
                return RANGE;
            case "fusion pair":
                return FUSION_PAIR;
            case "fusion promiscuous":
                return FUSION_PROMISCUOUS;
            case "Signatures":
                return SIGNATURE;
            case "MSI":
                return SIGNATURE_MSI;
            case "HRD":
                return SIGNATURE_HRD;
            case "MTL":
                return SIGNATURE_MTL;
            case "MTB":
                return SIGNATURE_MTB;
            default:
                LOGGER.warn("Unknown genomic event!");
                return UNKNOWN;
        }
    }
}
