package com.hartwig.hmftools.common.lims;

import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum LimsGermlineReportingLevel {
    REPORT_WITH_NOTIFICATION("report with notification"),
    REPORT_WITHOUT_NOTIFICATION("report without notification"),
    NO_REPORTING("no reporting");

    @NotNull
    private final String display;

    LimsGermlineReportingLevel(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }

    @NotNull
    static LimsGermlineReportingLevel fromLimsInputs(boolean limsSampleReportGermlineVariants, @NotNull String germlineReportingLevelString,
            @NotNull String sampleId, @Nullable LimsCohortConfig cohort) {
        if (limsSampleReportGermlineVariants && cohort != null) {
            // Cases "geen toevalsbevindingen: familie mag deze/wel niet opvragen" have been merged
            // into a single category "geen toevalsbevindingen" per feb 1st 2020
            if (cohort.reportGermline() && cohort.reportGermlineFlag()) {
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
                                "Cannot resolve germline reporting choice " + germlineReportingLevelString + " for sample " + sampleId);
                }
            } else if (cohort.reportGermline() && !cohort.reportGermlineFlag()) {
                return REPORT_WITHOUT_NOTIFICATION;
            }
        }
        return NO_REPORTING;
    }
}
