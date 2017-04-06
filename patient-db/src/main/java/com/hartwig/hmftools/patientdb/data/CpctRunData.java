package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import org.jetbrains.annotations.NotNull;

public class CpctRunData {
    @NotNull
    private final String folderName;
    @NotNull
    private final LocalDate uploadDate;
    @NotNull
    private final String referenceSampleId;
    @NotNull
    private final String tumorSampleId;
    @NotNull
    private final String patientId;

    public CpctRunData(@NotNull final String folderName, @NotNull final LocalDate uploadDate,
            @NotNull final String referenceSampleId, @NotNull final String tumorSampleId,
            @NotNull final String patientId) {
        this.uploadDate = uploadDate;
        this.patientId = patientId;
        this.referenceSampleId = referenceSampleId;
        this.tumorSampleId = tumorSampleId;
        this.folderName = folderName;
    }

    @NotNull
    public String patientId() {
        return patientId;
    }

    @NotNull
    public String folderName() {
        return folderName;
    }
}
