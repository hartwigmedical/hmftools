package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsGermlineReportingLevel {
    REPORT_WITH_NOTIFICATION,
    REPORT_WITHOUT_NOTIFICATION,
    NO_REPORTING;

    @NotNull
    static LimsGermlineReportingLevel fromLimsInputs(boolean reportGermlineVariants, @NotNull String germlineReportingLevelString,
            @NotNull String sampleId, @NotNull LimsCohort cohort) {
        if (reportGermlineVariants) {
            if (cohort == LimsCohort.WIDE) {
                // Cases "geen toevalsbevindingen: familie mag deze/wel niet opvragen" have been merged
                // into a single category "geen toevalsbevindingen" per feb 1st 2020
                switch (germlineReportingLevelString) {
                    case "1: Behandelbare toevalsbevindingen":
                    case "2: Alle toevalsbevindingen":
                        return REPORT_WITH_NOTIFICATION;
                    case "3: Geen toevalsbevindingen":
                    case "3: Geen toevalsbevindingen; familie mag deze wel opvragen":
                    case "4: Geen toevalsbevindingen; familie mag deze niet opvragen":
                        return REPORT_WITHOUT_NOTIFICATION;
                    default:
                        throw new IllegalStateException(
                                "Cannot resolve germline reporting choice for sample: " + sampleId + ": " + germlineReportingLevelString);
                }
            } else if (cohort == LimsCohort.COREDB) {
                switch (germlineReportingLevelString) {
                    case "1: Yes’":
                        return REPORT_WITH_NOTIFICATION;
                    case "2: No":
                        return REPORT_WITHOUT_NOTIFICATION;
                    default:
                        throw new IllegalStateException(
                                "Cannot resolve germline reporting choice for sample: " + sampleId + ": " + germlineReportingLevelString);
                }

            }
            return REPORT_WITHOUT_NOTIFICATION;
        }
        return NO_REPORTING;
    }
}
