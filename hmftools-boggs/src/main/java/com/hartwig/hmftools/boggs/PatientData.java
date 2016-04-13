package com.hartwig.hmftools.boggs;

import org.jetbrains.annotations.NotNull;

public class PatientData {
    @NotNull
    private final String pipelineRunDir;
    @NotNull
    private final String refSampleID;
    @NotNull
    private final String refExternalID;
    @NotNull
    private final String tumorSampleID;
    @NotNull
    private final String tumorExternalID;

    public PatientData(@NotNull String pipelineRunDir, @NotNull String refSampleID,
                       @NotNull String refExternalID, @NotNull String tumorSampleID, @NotNull String tumorExternalID) {
        this.pipelineRunDir = pipelineRunDir;
        this.refSampleID = refSampleID;
        this.refExternalID = refExternalID;
        this.tumorSampleID = tumorSampleID;
        this.tumorExternalID = tumorExternalID;
    }
}
