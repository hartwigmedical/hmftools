package com.hartwig.hmftools.serve.extraction.util;

public enum DriverInconsistencyMode {
    IGNORE(false),
    WARN_ONLY(true),
    FILTER(true);

    private final boolean active;

    DriverInconsistencyMode(final boolean active) {
        this.active = active;
    }

    public boolean isActive() {
        return active;
    }
}