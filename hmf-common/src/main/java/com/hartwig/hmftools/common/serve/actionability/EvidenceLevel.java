package com.hartwig.hmftools.common.serve.actionability;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum EvidenceLevel {
    A,
    B,
    C,
    D;

    @Nullable
    public static EvidenceLevel fromString(@NotNull String string) {
        switch (string) {
            case "A":
                return EvidenceLevel.A;
            case "B":
                return EvidenceLevel.B;
            case "C":
                return EvidenceLevel.C;
            case "D":
                return EvidenceLevel.D;
            case "NA":
            default:
                return null;
        }
    }
}
