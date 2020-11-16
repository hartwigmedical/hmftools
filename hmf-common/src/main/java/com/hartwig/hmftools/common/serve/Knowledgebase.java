package com.hartwig.hmftools.common.serve;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum Knowledgebase {
    CGI("CGI", false),
    CIVIC("CIViC", false),
    DOCM("DoCM", false),
    HARTWIG_COHORT("HartwigCohort", false),
    HARTWIG_CURATED("HartwigCurated", false),
    ICLUSION("iClusion", true),
    JAX("JAX", false),
    ONCOKB("OncoKB", false);

    @Nullable
    public static Knowledgebase fromDisplayString(@NotNull String display) {
        for (Knowledgebase knowledgebase : Knowledgebase.values()) {
            if (knowledgebase.display().equals(display)) {
                return knowledgebase;
            }
        }

        return null;
    }

    @NotNull
    private final String display;
    private final boolean trialSource;

    Knowledgebase(@NotNull final String display, final boolean trialSource) {
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
