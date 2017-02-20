package com.hartwig.hmftools.patientreporter.slicing;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

public class HMFSlicingAnnotation {

    @NotNull
    private final String transcriptID;
    private final int transcriptVersion;
    @NotNull
    private final String gene;

    HMFSlicingAnnotation(@NotNull final String transcriptID, final int transcriptVersion,
            @NotNull final String gene) {
        this.transcriptID = transcriptID;
        this.transcriptVersion = transcriptVersion;
        this.gene = gene;
    }

    @NotNull
    public String transcript() {
        return transcriptID + "." + transcriptVersion();
    }

    @NotNull
    public String transcriptID() {
        return transcriptID;
    }

    @VisibleForTesting
    int transcriptVersion() {
        return transcriptVersion;
    }

    @NotNull
    public String gene() {
        return gene;
    }
}
