package com.hartwig.hmftools.serve.actionability;

import com.hartwig.hmftools.common.serve.EvidenceDirection;
import com.hartwig.hmftools.common.serve.Source;

import org.jetbrains.annotations.NotNull;

public final class ActionableEventFactory {

    private ActionableEventFactory() {
    }

    @NotNull
    public static Source sourceFromFileValue(@NotNull String sourceFileValue) {
        Source source = Source.fromDisplayString(sourceFileValue);
        if (source == null) {
            throw new IllegalStateException("Cannot resolve source: " + sourceFileValue);
        }
        return source;
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
