package com.hartwig.hmftools.common.actionability;

import org.jetbrains.annotations.NotNull;

public enum ActionabilitySource {
    ICLUSION("iClusion", true),
    ONCOKB("OncoKb", false),
    CIVIC("CiViC", false),
    CGI("CGI", false);

    @NotNull
    private final String sourceName;

    private final boolean isTrialSource;

    ActionabilitySource(@NotNull final String sourceName, final boolean isTrialSource) {
        this.sourceName = sourceName;
        this.isTrialSource = isTrialSource;
    }

    @NotNull
    public String sourceName() {
        return sourceName;
    }

    public boolean isTrialSource() {
        return isTrialSource;
    }

    @NotNull
    public static ActionabilitySource fromString(@NotNull String source) {
        switch (source.toLowerCase()) {
            case "oncokb":
                return ONCOKB;
            case "cgi":
                return CGI;
            case "civic":
                return CIVIC;
            case "iclusion":
                return ICLUSION;
            default:
                throw new IllegalArgumentException("Unrecognized actionability source: " + source);
        }
    }
}
