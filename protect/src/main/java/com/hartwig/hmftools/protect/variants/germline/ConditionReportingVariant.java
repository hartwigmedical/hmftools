package com.hartwig.hmftools.protect.variants.germline;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum ConditionReportingVariant {
    BIALLELIC_ONLY("Biallelic"),
    ALL("Monoallelic"), // will reported both biallelic and monoallelic variants
    UNKNOWN("unknown");

    private static final Logger LOGGER = LogManager.getLogger(ConditionReportingVariant.class);

    @NotNull
    private final String display;

    ConditionReportingVariant(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public static ConditionReportingVariant fromConditionString(@NotNull String conditionString) {
        for (ConditionReportingVariant condition : ConditionReportingVariant.values()) {
            if (condition.display.equals(conditionString)) {
                return condition;
            }
        }

        LOGGER.warn("Unknown reporting condition: {}!", conditionString);
        return ConditionReportingVariant.UNKNOWN;
    }
}
