package com.hartwig.hmftools.serve.actionability;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum EvidenceDirection {
    RESPONSIVE("Responsive"),
    RESISTANT("Resistant");

    @Nullable
    public static EvidenceDirection fromDisplayString(@NotNull String display) {
        for (EvidenceDirection direction : EvidenceDirection.values()) {
            if (direction.display().equals(display)) {
                return direction;
            }
        }

        return null;
    }

    @NotNull
    private final String display;

    EvidenceDirection(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }
}
