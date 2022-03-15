package com.hartwig.hmftools.serve.extraction.catalog;

public enum DriverInconsistencyMode {
    IGNORE(false), //We want to report whole source
    WARN_ONLY(true), //We want report source when match with driver catalog
    FILTER(true); //We want report source when match with driver catalog + filtering

    private final boolean logging;

    DriverInconsistencyMode(final boolean logging) {
        this.logging = logging;
    }

    public boolean logging() {
        return logging;
    }
}