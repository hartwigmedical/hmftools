package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsGermlineReportingChoice {
    ALL,
    ACTIONABLE_ONLY,
    NONE_ALLOW_FAMILY,
    NONE,
    UNKNOWN;

    @NotNull
    static LimsGermlineReportingChoice fromLimsGermlineReportingChoiceString(@NotNull String germlineReportingChoiceString,
            @NotNull String sampleId) {
        LimsSampleType type = LimsSampleType.fromSampleId(sampleId);
        if (type == LimsSampleType.WIDE) {
            switch (germlineReportingChoiceString) {
                case "2: Alle toevalsbevindingen":
                    return ALL;
                case "1: Behandelbare toevalsbevindingen":
                    return ACTIONABLE_ONLY;
                case "3: Geen toevalsbevindingen; familie mag deze wel opvragen":
                    return NONE_ALLOW_FAMILY;
                case "3: Geen toevalsbevindingen":
                    return NONE;
                case "4: Geen toevalsbevindingen; familie mag deze niet opvragen":
                    return NONE;
                default:
                    throw new IllegalStateException(
                            "Cannot resolve germline reporting choice for sample: " + sampleId + ": " + germlineReportingChoiceString);
            }
        }

        return UNKNOWN;
    }
}
