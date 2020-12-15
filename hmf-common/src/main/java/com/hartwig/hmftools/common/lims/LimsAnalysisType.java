package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsAnalysisType {
    SOMATIC_T,
    SOMATIC_R;

    @NotNull
    public static LimsAnalysisType extractAnalysisType(@NotNull String analysisType) {
        switch (analysisType) {
            case "Somatic_T":
                return SOMATIC_T;
            case "Somatic_R":
                return SOMATIC_R;
            default:
                throw new IllegalStateException("Cannot resolve analysis type '{}' " + analysisType);
        }

    }
}
