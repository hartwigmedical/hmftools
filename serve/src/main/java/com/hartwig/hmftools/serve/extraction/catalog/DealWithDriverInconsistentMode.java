package com.hartwig.hmftools.serve.extraction.catalog;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class DealWithDriverInconsistentMode {

    private static final Logger LOGGER = LogManager.getLogger(DealWithDriverInconsistentMode.class);

    public static boolean filterOnInconsistenties(@NotNull DealWithDriverInconsistentModeAnnotation annotation) {
        if (annotation.equals(DealWithDriverInconsistentModeAnnotation.FILTER)) {
            return true;
        } else if (annotation.equals(DealWithDriverInconsistentModeAnnotation.IGNORE)) {
            return true;
        } else if (annotation.equals(DealWithDriverInconsistentModeAnnotation.WARN_ONLY)) {
            return false;
        } else {
            return false;
        }
    }

}
