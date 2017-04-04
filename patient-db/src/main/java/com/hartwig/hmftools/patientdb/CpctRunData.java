package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;

import org.jetbrains.annotations.NotNull;

class CpctRunData {
    @NotNull
    private final LocalDate uploadDate;
    @NotNull
    private final String referenceSampleId;
    @NotNull
    private final String tumorSampleId;
    @NotNull
    private final String patientId;

    CpctRunData(@NotNull final LocalDate uploadDate, @NotNull final String referenceSampleId,
            @NotNull final String tumorSampleId, @NotNull final String patientId) {
        this.uploadDate = uploadDate;
        this.patientId = patientId;
        this.referenceSampleId = referenceSampleId;
        this.tumorSampleId = tumorSampleId;
    }

    @Override
    public String toString() {
        final StringBuffer bf = new StringBuffer();
        bf.append(uploadDate).append("-").append(patientId).append("-").append(referenceSampleId).append("-").append(
                tumorSampleId).append("\n");
        return bf.toString();
    }

    @NotNull
    String patientId() {
        return patientId;
    }
}
