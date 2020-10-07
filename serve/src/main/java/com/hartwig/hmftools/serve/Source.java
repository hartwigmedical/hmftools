package com.hartwig.hmftools.serve;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum Source {
    CIVIC("CIViC"),
    ONCOKB("OncoKB"),
    CGI("CGI"),
    JAX("JAX"),
    ICLUSION("iClusion");

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
