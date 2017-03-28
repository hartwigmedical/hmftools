package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;

import org.jetbrains.annotations.Nullable;

class TreatmentData {
    //    private final String type;
    private final String treatmentName;
    private final LocalDate startDate;
    private final LocalDate endDate;
    private final String earlyResponse;

    TreatmentData(@Nullable LocalDate startDate, @Nullable LocalDate endDate, @Nullable String treatmentName,
            @Nullable String earlyResponse) {
        this.startDate = startDate;
        this.endDate = endDate;
        this.treatmentName = treatmentName;
        this.earlyResponse = earlyResponse;
    }

    @Override
    public String toString() {
        final StringBuffer bf = new StringBuffer();
        bf.append("(").append(startDate).append("->").append(endDate).append("): ").append(treatmentName).append(
                " - ").append(earlyResponse).append("\n");
        return bf.toString();
    }

    String treatmentName() {
        return treatmentName;
    }

    LocalDate startDate() {
        return startDate;
    }

    LocalDate endDate() {
        return endDate;
    }

    String earlyResponse() {
        return earlyResponse;
    }
}
