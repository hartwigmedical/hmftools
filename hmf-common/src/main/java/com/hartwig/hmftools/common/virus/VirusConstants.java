package com.hartwig.hmftools.common.virus;

import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public enum VirusConstants {
    MCV("MCV", true),
    EBV("EBV", true),
    HPV("HPV", true),
    HBV("HBV", false),
    HHV8("HHV-8", false);

    @NotNull
    private final String virusName;
    private final boolean reportVirusOnSummary;

    VirusConstants(@NotNull final String virusName, final boolean reportVirusOnSummary) {
        this.virusName = virusName;
        this.reportVirusOnSummary = reportVirusOnSummary;
    }

    public boolean reportVirusOnSummary() {
        return reportVirusOnSummary;
    }

    @NotNull
    public static VirusConstants fromVirusName(@NotNull String virusName) {
        switch (virusName) {
            case "MCV":
                return MCV;
            case "EBV":
                return EBV;
            case "HPV":
                return HPV;
            case "HBV":
                return HBV;
            case "HHV-8":
                return HHV8;
            default:
                throw new IllegalStateException("Cannot resolve virus name: " + virusName);
        }
    }

    @NotNull
    public static List<String> allViruses() {
        List<String> virusSummary = Lists.newArrayList();
        for (VirusConstants virus : VirusConstants.values()) {
            if (virus.reportVirusOnSummary()) {
                virusSummary.add(virus.virusName);
            }
        }

        return virusSummary;
    }
}
