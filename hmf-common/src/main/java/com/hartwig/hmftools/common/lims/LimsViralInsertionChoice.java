package com.hartwig.hmftools.common.lims;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum LimsViralInsertionChoice {
    REPORT_VIRAL_INSERION,
    NO_REPORT_VIRAL_INSERTIONS;

    private static final Logger LOGGER = LogManager.getLogger(LimsViralInsertionChoice.class);

    @NotNull
    static LimsViralInsertionChoice fromLimsViralInsertionsReportingChoiceString(boolean viralInsertionsReportingChoiceString,
            @NotNull String sampleId) {
        LimsStudy study = LimsStudy.fromSampleId(sampleId);

        if (viralInsertionsReportingChoiceString) {
            if (study == LimsStudy.DRUP || study == LimsStudy.CPCT) {
                LOGGER.warn("Consent of viral insertions are true, but must be false for CPCT/DRUP");
            }
            return REPORT_VIRAL_INSERION;
        } else {
            if (study == LimsStudy.CORE || study == LimsStudy.WIDE) {
                LOGGER.warn("Consent of viral insertions are false, but must be true for WIDE/CORE");
            }
            return NO_REPORT_VIRAL_INSERTIONS;
        }
    }
}
