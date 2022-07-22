package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;

public enum LimsAnalysisType {
    SOMATIC_T,
    SOMATIC_R,
    TARGETED_TUMOR_ONLY;

    @NotNull
    public static LimsAnalysisType extractAnalysisType(@NotNull String analysisType) {
        switch (analysisType) {
            case "Somatic_T":
                return SOMATIC_T;
            case "Somatic_R":
                return SOMATIC_R;
            case "Targeted_Tumor_Only":
                return TARGETED_TUMOR_ONLY;
            default:
                throw new IllegalStateException("Cannot resolve analysis type '{}' " + analysisType);
        }
    }
}
