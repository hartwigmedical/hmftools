package com.hartwig.hmftools.common.lims;

import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public enum LimsGermlineReportingLevel {
    REPORT_WITH_NOTIFICATION,
    REPORT_WITHOUT_NOTIFICATION,
    NO_REPORTING;

    private static final Logger LOGGER = LogManager.getLogger(LimsGermlineReportingLevel.class);

    @NotNull
    static LimsGermlineReportingLevel fromLimsInputs(boolean reportGermlineVariants, @NotNull String germlineReportingLevelString,
            @NotNull String sampleId, @NotNull LimsCohortConfig cohort) {
        if (reportGermlineVariants && !cohort.cohortId().equals(Strings.EMPTY)) {
            // Cases "geen toevalsbevindingen: familie mag deze/wel niet opvragen" have been merged
            // into a single category "geen toevalsbevindingen" per feb 1st 2020
            if (cohort.reportGermline() || cohort.reportGermlineFlag()) {
                switch (germlineReportingLevelString) {
                    case "1: Behandelbare toevalsbevindingen":
                    case "2: Alle toevalsbevindingen":
                    case "1: Yes":
                        return REPORT_WITH_NOTIFICATION;
                    case "3: Geen toevalsbevindingen":
                    case "3: Geen toevalsbevindingen; familie mag deze wel opvragen":
                    case "4: Geen toevalsbevindingen; familie mag deze niet opvragen":
                    case "2: No":
                        return REPORT_WITHOUT_NOTIFICATION;
                    default:
                        throw new IllegalStateException(
                                "Cannot resolve germline reporting choice " + germlineReportingLevelString + " for sample " + sampleId );
                }
            } else if (cohort.reportGermline() || !cohort.reportGermlineFlag()) {

                return REPORT_WITHOUT_NOTIFICATION;
            }
        }
        return NO_REPORTING;
    }
}
