package com.hartwig.hmftools.serve.actionability;

import com.hartwig.hmftools.common.serve.EvidenceDirection;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public final class ActionableEventFactory {

    private ActionableEventFactory() {
    }

    @NotNull
    public static Knowledgebase sourceFromFileValue(@NotNull String sourceFileValue) {
        Knowledgebase source = Knowledgebase.fromDisplayString(sourceFileValue);
        if (source == null) {
            throw new IllegalStateException("Cannot resolve source knowledgebase: " + sourceFileValue);
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
