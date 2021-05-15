package com.hartwig.hmftools.patientreporter.germline;

import org.jetbrains.annotations.NotNull;

public enum GermlineCondition {
    NEVER,
    ALWAYS,
    ONLY_GERMLINE_HOM,
    ONLY_SPECIFIC_VARIANT;

    @NotNull
    static GermlineCondition toGermlineCondition(@NotNull String germlineConditionInput) {
        for (GermlineCondition condition : GermlineCondition.values()) {
            if (germlineConditionInput.equals(condition.toString())) {
                return condition;
            }
        }

        throw new IllegalStateException("Cannot resolve germline condition: '{}' " + germlineConditionInput);
    }
}
