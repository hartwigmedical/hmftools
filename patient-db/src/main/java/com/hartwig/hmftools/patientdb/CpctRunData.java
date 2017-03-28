package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;

import org.jetbrains.annotations.NotNull;

class CpctRunData {
    private final LocalDate uploadDate;
    private final String referenceSampleId;
    private final String tumorSampleId;
    private final String patientId;

    CpctRunData(@NotNull LocalDate uploadDate, @NotNull String referenceSampleId, @NotNull String tumorSampleId,
            @NotNull String patientId) {
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
