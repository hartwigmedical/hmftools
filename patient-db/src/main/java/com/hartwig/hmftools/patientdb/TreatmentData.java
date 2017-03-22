package com.hartwig.hmftools.patientdb;

import java.util.Date;

import org.jetbrains.annotations.Nullable;

public class TreatmentData {
    //    private final String type;
    private final String treatmentName;
    private final Date startDate;
    private final Date endDate;
    private final String earlyResponse;

    TreatmentData(@Nullable Date startDate, @Nullable Date endDate, @Nullable String treatmentName,
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
}
