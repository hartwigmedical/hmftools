package com.hartwig.hmftools.serve.extraction.catalog;

import org.jetbrains.annotations.NotNull;

public class DealWithDriverInconsistentMode {

    public static boolean filterOnInconsistenties(@NotNull DealWithDriverInconsistentModeAnnotation annotation) {
        switch (annotation) {
            case FILTER:
            case WARN_ONLY:
                return true;
            case IGNORE:
                return false;
            default:
                return false;
        }
    }
}