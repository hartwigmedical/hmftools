package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsGermlineReportingChoice {
    ALL,
    ACTIONABLE_ONLY,
    NONE,
    UNKNOWN;

    @NotNull
    static LimsGermlineReportingChoice fromLimsGermlineReportingChoiceString(@NotNull String germlineReportingChoiceString,
            @NotNull String sampleId) {
        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);
        // Cases "geen toevalsbevindingen: familie mag deze/wel niet opvragen" have been merged
        // into a single category "geen toevalsbevindingen" per feb 1st 2020
        if (type == LimsSampleType.WIDE) {
            switch (germlineReportingChoiceString) {
                case "2: Alle toevalsbevindingen":
                    return ALL;
                case "1: Behandelbare toevalsbevindingen":
                    return ACTIONABLE_ONLY;
                case "3: Geen toevalsbevindingen":
                case "3: Geen toevalsbevindingen; familie mag deze wel opvragen":
                case "4: Geen toevalsbevindingen; familie mag deze niet opvragen":
                    return NONE;
                default:
                    throw new IllegalStateException(
                            "Cannot resolve germline reporting choice for sample: " + sampleId + ": " + germlineReportingChoiceString);
            }
        } else if (type == LimsSampleType.CORE) {
            return NONE;
        }

        return UNKNOWN;
    }
}
