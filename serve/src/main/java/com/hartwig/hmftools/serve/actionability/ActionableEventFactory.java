package com.hartwig.hmftools.serve.actionability;

import org.jetbrains.annotations.NotNull;

public final class ActionableEventFactory {

    private ActionableEventFactory() {
    }

    @NotNull
    public static EvidenceDirection directionFromFileValue(@NotNull String directionFileValue) {
        EvidenceDirection direction = EvidenceDirection.fromDisplayString(directionFileValue);
        if (direction == null) {
            throw new IllegalStateException("Cannot resolve direction: " + directionFileValue);
        }
        return direction;
    }
}
