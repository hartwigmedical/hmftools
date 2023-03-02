package com.hartwig.hmftools.datamodel.genotype;

import org.jetbrains.annotations.NotNull;

public enum GenotypeStatus {
    HOM_REF("HOM"),
    HET("HET"),
    HOM_ALT("HOM"),
    UNKNOWN("UNKNOWN");

    @NotNull
    private final String simpleDisplay;

    GenotypeStatus(@NotNull final String simpleDisplay) {
        this.simpleDisplay = simpleDisplay;
    }

    @NotNull
    public String simplifiedDisplay() {
        return simpleDisplay;
    }
}
