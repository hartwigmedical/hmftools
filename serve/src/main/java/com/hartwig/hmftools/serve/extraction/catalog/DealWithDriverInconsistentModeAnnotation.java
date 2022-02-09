package com.hartwig.hmftools.serve.extraction.catalog;

import org.jetbrains.annotations.NotNull;

public enum DealWithDriverInconsistentModeAnnotation {
    IGNORE, //We want to report whole source
    WARN_ONLY, //We want report source when match with driver catalog
    FILTER, //We want report source when match with driver catalog + filtering
    UNKNOWN;

    @NotNull
    public static DealWithDriverInconsistentModeAnnotation extractDealWithDriverInconsistentMode(
            @NotNull String dealWithDriverInconsistentMode) {
        if (dealWithDriverInconsistentMode.equals("ignore")) {
            return IGNORE;
        } else if (dealWithDriverInconsistentMode.equals("warn_only")) {
            return WARN_ONLY;
        } else if (dealWithDriverInconsistentMode.equals("filter")) {
            return FILTER;
        } else {
            return UNKNOWN;
        }
    }
}
