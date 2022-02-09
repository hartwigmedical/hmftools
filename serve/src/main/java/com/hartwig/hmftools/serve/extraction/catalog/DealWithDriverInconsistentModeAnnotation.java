package com.hartwig.hmftools.serve.extraction.catalog;

import org.jetbrains.annotations.NotNull;

public enum DealWithDriverInconsistentModeAnnotation {
    IGNORE, //We want not to deal with any driver inconsistency!
    WARN_ONLY, //We want to deal with any driver inconsistency!
    FILTER, //We want to deal with any driver inconsistency! + remove
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
