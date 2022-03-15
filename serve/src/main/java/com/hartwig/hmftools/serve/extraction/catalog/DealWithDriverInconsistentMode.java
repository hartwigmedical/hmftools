package com.hartwig.hmftools.serve.extraction.catalog;

import org.jetbrains.annotations.NotNull;

public final class DealWithDriverInconsistentMode {

    public static boolean filterOnInconsistencies(@NotNull DriverInconsistencyMode mode) {
        switch (mode) {
            case FILTER:
            case WARN_ONLY:
                return true;
            default:
                return false;
        }
    }
}