package com.hartwig.hmftools.patientreporter.variants;

import org.jetbrains.annotations.NotNull;

public enum ClonalInterpretation {
    LIKELY("Likely"),
    UNCERTAIN("Uncertain"),
    UNLIKELY("Unlikely");

    @NotNull
    private final String text;

    ClonalInterpretation(@NotNull final String text) {
        this.text = text;
    }

    @NotNull
    public String text() {
        return text;
    }

    @NotNull
    public static ClonalInterpretation interpret(double clonalLikelihood) {
        if (clonalLikelihood > 0.8) {
            return LIKELY ;
        } else if (clonalLikelihood > 0.2) {
            return UNCERTAIN;
        } else {
            return UNLIKELY;
        }
    }
}
