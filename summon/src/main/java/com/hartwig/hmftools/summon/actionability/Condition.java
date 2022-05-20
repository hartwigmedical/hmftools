package com.hartwig.hmftools.summon.actionability;

import org.jetbrains.annotations.NotNull;

public enum Condition {
    ONLY_HIGH,
    ALWAYS,
    ALWAYS_NO_ACTIONABLE,
    OTHER;

    @NotNull
    static Condition toCondition(@NotNull String conditionInput) {
        for (Condition condition : Condition.values()) {
            if (conditionInput.equals(condition.toString())) {
                return condition;
            }
        }

        throw new IllegalStateException("Cannot resolve condition: '{}' " + conditionInput);
    }
}
