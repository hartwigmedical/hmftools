package com.hartwig.hmftools.healthchecker.runners.checks;

import org.jetbrains.annotations.NotNull;

public enum SomaticCheck {
    COUNT_TOTAL("SOMATIC_%s_COUNT"),
    DBSNP_COUNT("SOMATIC_%s_DBSNP_COUNT"),
    AF_MEDIAN("SOMATIC_AF_MEDIAN"),
    AF_LOWER_SD("SOMATIC_AF_LOWER_SD"),
    AF_UPPER_SD("SOMATIC_AF_UPPER_SD"),
    PROPORTION_CHECK("SOMATIC_%s_PROPORTION_VARIANTS");

    @NotNull
    private final String checkPattern;

    SomaticCheck(@NotNull final String label) {
        this.checkPattern = label;
    }

    @NotNull
    public String checkName() {
        return checkPattern;
    }

    @NotNull
    public String checkName(@NotNull final String input1) {
        return String.format(checkPattern, input1.toUpperCase());
    }

    @NotNull
    public String checkName(@NotNull final String input1, @NotNull final String input2) {
        return String.format(checkPattern, input1.toUpperCase(), input2.toUpperCase());
    }
}
