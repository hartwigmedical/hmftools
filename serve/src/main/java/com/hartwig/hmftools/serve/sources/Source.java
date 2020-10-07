package com.hartwig.hmftools.serve.sources;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum Source {
    CGI("CGI"),
    CIVIC("CIViC"),
    DOCM("DoCM"),
    HARTWIG_COHORT("HartwigCohort"),
    HARTWIG_CURATED("HartwigCurated"),
    ICLUSION("iClusion"),
    JAX("JAX"),
    ONCOKB("OncoKB");

    @Nullable
    public static Source fromDisplayString(@NotNull String display) {
        for (Source source : Source.values()) {
            if (source.display().equals(display)) {
                return source;
            }
        }

        return null;
    }

    @NotNull
    private final String display;

    Source(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }
}
