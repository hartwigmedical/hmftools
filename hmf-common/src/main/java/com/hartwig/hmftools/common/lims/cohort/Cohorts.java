package com.hartwig.hmftools.common.lims.cohort;

import org.jetbrains.annotations.NotNull;

public enum Cohorts {
    CORE("CORE"),
    WIDE("WIDE");

    @NotNull
    private final String display;

    Cohorts(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }
}
