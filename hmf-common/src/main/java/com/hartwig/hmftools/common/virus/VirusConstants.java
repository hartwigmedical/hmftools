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

    private final boolean reportedVirusOnSummary;

    VirusConstants(@NotNull final String virusName, final boolean reportedSummary) {
        this.virusName = virusName;
        this.reportedVirusOnSummary = reportedSummary;
    }

    public boolean reportedVirusOnSummary() {
        return reportedVirusOnSummary;
    }

    @NotNull
    public static VirusConstants virusName(@NotNull String virusName) {

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
                throw new IllegalStateException("Cannot resolve firus name: " + virusName);
        }
    }

    @NotNull
    public static List<String> allVirussen() {
        List<String> virusSummary = Lists.newArrayList();
        for (VirusConstants virussen : VirusConstants.values()) {
            if (virussen.reportedVirusOnSummary) {
                virusSummary.add(virussen.virusName);
            }
        }

        return virusSummary;
    }
}
