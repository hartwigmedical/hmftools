package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsGermlineReportingChoice {
    REPORT_WITH_NOTIFICATION,
    REPORT_WITHOUT_NOTIFICATION,
    NO_REPORTING;

    @NotNull
    static LimsGermlineReportingChoice fromLimsGermlineReportingChoiceString(@NotNull String germlineReportingChoiceString,
            @NotNull String sampleId) {
        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);
        // Cases "geen toevalsbevindingen: familie mag deze/wel niet opvragen" have been merged
        // into a single category "geen toevalsbevindingen" per feb 1st 2020
        if (type == LimsSampleType.WIDE) {
            switch (germlineReportingChoiceString) {
                case "1: Behandelbare toevalsbevindingen":
                case "2: Alle toevalsbevindingen":
                    return REPORT_WITH_NOTIFICATION;
                case "3: Geen toevalsbevindingen":
                case "3: Geen toevalsbevindingen; familie mag deze wel opvragen":
                case "4: Geen toevalsbevindingen; familie mag deze niet opvragen":
                    return REPORT_WITHOUT_NOTIFICATION;
                default:
                    throw new IllegalStateException(
                            "Cannot resolve germline reporting choice for sample: " + sampleId + ": " + germlineReportingChoiceString);
            }
        } else if (type == LimsSampleType.CORE) {
            return REPORT_WITHOUT_NOTIFICATION;
        }

        return NO_REPORTING;
    }
}
