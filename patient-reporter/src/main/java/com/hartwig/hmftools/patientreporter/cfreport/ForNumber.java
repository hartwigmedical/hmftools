package com.hartwig.hmftools.patientreporter.cfreport;

import org.jetbrains.annotations.NotNull;

public enum ForNumber {
    FOR_082("082"),
    FOR_083("083"),
    FOR_100("100"),
    FOR_102("102"),
    FOR_080("080"),
    FOR_103("XXX"),
    FOR_UNDEFINED("Undefined");

    @NotNull
    public final String forNumber;

    ForNumber(@NotNull final String forNumber) {
        this.forNumber = forNumber;
    }

    @NotNull
    public String display() {
        return forNumber;
    }

}
