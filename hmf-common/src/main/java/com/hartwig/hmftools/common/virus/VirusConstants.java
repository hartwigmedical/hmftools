package com.hartwig.hmftools.common.virus;

import org.jetbrains.annotations.NotNull;

public enum VirusConstants {
    MCV("MCV"),
    EBV("EBV"),
    HPV("HPV"),
    HBV("HBV"),
    HHV8("HHV-8");

    @NotNull
    private final String virusName;

    VirusConstants(@NotNull final String virusName) {
        this.virusName = virusName;
    }

    @NotNull
    public static VirusConstants fromVirusName(@NotNull String virusName) {
        switch (virusName) {
            case "MCV":
                return MCV;
            case "EBV":
                return EBV;
            case "HPV":
                return HPV;
            case "HBV":
                return HBV;
            case "HHV-8":
                return HHV8;
            default:
                throw new IllegalStateException("Cannot resolve virus name: " + virusName);
        }
    }
}