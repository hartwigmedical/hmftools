package com.hartwig.hmftools.common.region.hmfslicer;

import org.jetbrains.annotations.NotNull;

public enum Strand {
    FORWARD,
    REVERSE;

    @NotNull
    public static Strand valueOf(int direction) {
        switch (direction) {
            case 1:
                return Strand.FORWARD;
            case -1:
                return Strand.REVERSE;
        }

        throw new IllegalArgumentException("Invalid direction " + direction);
    }
}
