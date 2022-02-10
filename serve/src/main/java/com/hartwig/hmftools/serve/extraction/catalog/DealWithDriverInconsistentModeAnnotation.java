package com.hartwig.hmftools.serve.extraction.catalog;

import org.jetbrains.annotations.NotNull;

public enum DealWithDriverInconsistentModeAnnotation {
    IGNORE(true), //We want to report whole source
    WARN_ONLY(true), //We want report source when match with driver catalog
    FILTER(true); //We want report source when match with driver catalog + filtering

    private final boolean logging;

    DealWithDriverInconsistentModeAnnotation(final boolean logging) {
        this.logging = logging;
    }

    public boolean logging() {
        return logging;
    }

    @NotNull
    public static DealWithDriverInconsistentModeAnnotation extractDealWithDriverInconsistentMode(
            @NotNull String dealWithDriverInconsistentMode) {

        switch (dealWithDriverInconsistentMode) {
            case "ignore":
                return IGNORE;
            case "warn_only":
                return WARN_ONLY;
            case "filter":
                return FILTER;
            default:
                throw new IllegalStateException(
                        "Cannot resolve driver inconsistent mode of " + dealWithDriverInconsistentMode);
        }
    }
}