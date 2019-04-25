package com.hartwig.hmftools.common.lims;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum LimsGermlineFindingsChoice {
    ALL_ACTIONABLE,
    ALL,
    NONE_FAMILY,
    NONE,
    UNKNOWN;

    private static final Logger LOGGER = LogManager.getLogger(LimsGermlineFindingsChoice.class);

    @NotNull
    public static LimsGermlineFindingsChoice extractChoiceInformedConsent(@NotNull String germlineOptions, @NotNull String sampleId) {
        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);
        if (type == LimsSampleType.WIDE) {
            switch (germlineOptions) {
                case "1: Behandelbare toevalsbevindingen":
                    LOGGER.info("Patient has given consent for: 1: Behandelbare toevalsbevindingen");
                    return ALL_ACTIONABLE;
                case "2: Alle toevalsbevindingen":
                    LOGGER.info("Patient has given consent for: 2: Alle toevalsbevindingen");
                    return ALL;
                case "3: Geen toevalsbevindingen, familie mag deze wel opvragen":
                    LOGGER.info("Patient has given consent for: 3: Geen toevalsbevindingen, familie mag deze wel opvragen");
                    return NONE_FAMILY;
                case "4: Geen toevalsbevindingen, familie mag deze niet opvragen":
                    LOGGER.info("Patient has given consent for: 4: Geen toevalsbevindingen, familie mag deze niet opvragen");
                    return NONE;
                default:
                    throw new IllegalStateException("Cannot resolve germline choice for sampleId: " + sampleId + ": " + germlineOptions);
            }
        } else {
            LOGGER.info("SampleId is not a sample of WIDE study!");
            return UNKNOWN;
        }
    }
}
