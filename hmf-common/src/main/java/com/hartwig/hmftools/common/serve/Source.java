package com.hartwig.hmftools.common.serve;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum Source {
    CGI("CGI", false),
    CIVIC("CIViC", false),
    DOCM("DoCM", false),
    HARTWIG_COHORT("HartwigCohort", false),
    HARTWIG_CURATED("HartwigCurated", false),
    ICLUSION("iClusion", true),
    JAX("JAX", false),
    ONCOKB("OncoKB", false);

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
    private final boolean trialSource;

    Source(@NotNull final String display, final boolean trialSource) {
        this.display = display;
        this.trialSource = trialSource;
    }

    @NotNull
    public String display() {
        return display;
    }

    public boolean isTrialSource() {
        return trialSource;
    }
}
