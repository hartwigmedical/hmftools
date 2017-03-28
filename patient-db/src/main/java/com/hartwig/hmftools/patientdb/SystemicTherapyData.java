package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;

import org.jetbrains.annotations.Nullable;

public class SystemicTherapyData {
    private final LocalDate startDate;
    private final LocalDate endDate;
    private final String type;
    private final String treatment;
    private final String bestResponse;

    SystemicTherapyData(@Nullable LocalDate startDate, @Nullable LocalDate endDate, @Nullable String type,
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

    LocalDate endDate() {
        return endDate;
    }

    LocalDate startDate() {
        return startDate;
    }

    String bestResponse() {
        return bestResponse;
    }

    String treatment() {
        return treatment;
    }

    String type() {
        return type;
    }
}
