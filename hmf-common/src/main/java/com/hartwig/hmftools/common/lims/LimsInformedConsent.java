package com.hartwig.hmftools.common.lims;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum LimsInformedConsent {
    ALL_ACTIONABLE,
    ALL,
    NONE_FAMILY,
    NONE,
    UNKNOWN;

    private static final Logger LOGGER = LogManager.getLogger(LimsInformedConsent.class);

    @NotNull
    public static LimsInformedConsent extractChoiceInformedConsent(@NotNull String germlineOptions, @NotNull String sampleId) {
        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);
        if (type == LimsSampleType.WIDE) {
            switch (germlineOptions) {
                case "1: Behandelbare toevalsbevindingen":
                    LOGGER.info("Patient has given informed consent for: 1: Behandelbare toevalsbevindingen");
                    return ALL_ACTIONABLE;
                case "2: Alle toevalsbevindingen":
                    LOGGER.info("Patient has given informed consent for: 2: Alle toevalsbevindingen");
                    return ALL;
                case "3: Geen toevalsbevindingen, familie mag deze wel opvragen":
                    LOGGER.info("Patient has given informed consent for: 3: Geen toevalsbevindingen, familie mag deze wel opvragen");
                    return NONE_FAMILY;
                case "4: Geen toevalsbevindingen, familie mag deze niet opvragen":
                    LOGGER.info("Patient has given informed consent for: 4: Geen toevalsbevindingen, familie mag deze niet opvragen");
                    return NONE;
                default:
                    throw new IllegalStateException("Cannot resolve choice germline informed consent for sampleId: " + sampleId);
            }
        } else {
            LOGGER.info("SampleId is not a sample of WIDE study!");
            return UNKNOWN;
        }
    }
}
