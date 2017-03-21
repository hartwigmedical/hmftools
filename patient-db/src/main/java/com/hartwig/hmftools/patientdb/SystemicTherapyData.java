package com.hartwig.hmftools.patientdb;

import java.util.Date;

import org.jetbrains.annotations.Nullable;

public class SystemicTherapyData {
    private final Date startDate;
    private final Date endDate;
    private final String type;
    private final String treatment;
    private final String bestResponse;

    SystemicTherapyData(@Nullable Date startDate, @Nullable Date endDate, @Nullable String type,
            @Nullable String treatment, @Nullable String bestResponse) {
        this.startDate = startDate;
        this.endDate = endDate;
        this.type = type;
        this.treatment = treatment;
        this.bestResponse = bestResponse;
    }

    @Override
    public String toString() {
        final StringBuffer bf = new StringBuffer();
        bf.append("(").append(startDate).append("->").append(endDate).append("): ").append(type).append(" - ").append(
                treatment).append(" - ").append(bestResponse).append("\n");
        return bf.toString();
    }
}
