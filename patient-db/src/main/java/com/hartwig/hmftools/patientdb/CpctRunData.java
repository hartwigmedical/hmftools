package com.hartwig.hmftools.patientdb;

import java.util.Date;

import org.jetbrains.annotations.NotNull;

public class CpctRunData {
    private final Date uploadDate;
    private final String referenceSampleId;
    private final String tumorSampleId;
    private final String patientId;

    CpctRunData(@NotNull Date uploadDate, @NotNull String referenceSampleId, @NotNull String tumorSampleId,
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
