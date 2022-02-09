package com.hartwig.hmftools.serve.extraction.catalog;

import org.jetbrains.annotations.NotNull;

public class DealWithDriverInconsistentMode {

    public static boolean filterOnInconsistenties(@NotNull DealWithDriverInconsistentModeAnnotation annotation) {
        if (annotation.equals(DealWithDriverInconsistentModeAnnotation.FILTER)) {
            return false;
        } else if (annotation.equals(DealWithDriverInconsistentModeAnnotation.IGNORE)) {
            return true;
        } else if (annotation.equals(DealWithDriverInconsistentModeAnnotation.WARN_ONLY)) {
            return false;
        } else {
            return false;
        }
    }
}