package com.hartwig.hmftools.common.chromosome;

import org.jetbrains.annotations.NotNull;

public final class Chromosomes {

    private Chromosomes() {
    }

    @Deprecated
    public static int asInt(@NotNull final String chromosome) {
        switch (chromosome) {
            case "X": return 23;
            case "Y": return 24;
            case "MT": return 25;
        }
        return Integer.valueOf(chromosome);
    }

}
