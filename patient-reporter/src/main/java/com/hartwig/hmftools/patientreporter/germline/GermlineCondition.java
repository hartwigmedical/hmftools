package com.hartwig.hmftools.patientreporter.germline;


import org.jetbrains.annotations.NotNull;

public enum GermlineCondition {
    NEVER,
    ALWAYS,
    ONLY_GERMLINE_HOM,
    ONLY_SPECIFIC_VARIANT;

    @NotNull
    public static GermlineCondition extractGermlineCondition(@NotNull String germlineConditionInput) {
        switch (germlineConditionInput) {
            case "ALWAYS":
                return ALWAYS;
            case "NEVER":
                return NEVER;
            case "ONLY_GERMLINE_HOM":
                return ONLY_GERMLINE_HOM;
            case "ONLY_SPECIFIC_VARIANT":
                return ONLY_SPECIFIC_VARIANT;
            default:
                throw new IllegalStateException("Cannot resolve germline condition notify '{}' " + germlineConditionInput);
        }
    }
}
