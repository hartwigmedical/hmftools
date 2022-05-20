package com.hartwig.hmftools.patientreporter;

import org.jetbrains.annotations.NotNull;

public enum QsFormNumber {
    FOR_082("HMF-FOR-082"),
    FOR_083("HMF-FOR-083"),
    FOR_100("HMF-FOR-100"),
    FOR_102("HMF-FOR-102"),
    FOR_080("HMF-FOR-080"),
    FOR_209("HMF-FOR-209"),
    FOR_344("HMF-FOR-344"),
    FOR_345("HMF-FOR-345");

    @NotNull
    private final String display;

    QsFormNumber(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }

}
