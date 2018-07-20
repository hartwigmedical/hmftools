package com.hartwig.hmftools.breakpointinspector.datamodel;

import org.jetbrains.annotations.NotNull;

public enum HMFStructuralVariantType {
    INS(1, -1, "INNIE", false),
    DEL(1, -1, "INNIE", false),
    INV3(1, 1, "TANDEM_RIGHT", true),
    INV5(-1, -1, "TANDEM_LEFT", true),
    DUP(-1, 1, "OUTIE", false);

    private final int orientationBP1;
    private final int orientationBP2;
    @NotNull
    private final String orientationString;
    private final boolean inversion;

    HMFStructuralVariantType(final int orientationBP1, final int orientationBP2, @NotNull final String orientationString,
            final boolean inversion) {
        this.orientationBP1 = orientationBP1;
        this.orientationBP2 = orientationBP2;
        this.orientationString = orientationString;
        this.inversion = inversion;
    }

    public int orientationBP1() {
        return orientationBP1;
    }

    public int orientationBP2() {
        return orientationBP2;
    }

    @NotNull
    public String orientationString() {
        return orientationString;
    }

    public boolean isInversion() {
        return inversion;
    }
}
