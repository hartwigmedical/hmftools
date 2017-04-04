package com.hartwig.hmftools.patientdb;

import java.time.LocalDate;

import org.jetbrains.annotations.Nullable;

public class SystemicTherapyData {
    @Nullable
    private final LocalDate startDate;
    @Nullable
    private final LocalDate endDate;
    @Nullable
    private final String type;
    @Nullable
    private final String treatment;
    @Nullable
    private final String bestResponse;

    SystemicTherapyData(@Nullable final LocalDate startDate, @Nullable final LocalDate endDate,
            @Nullable final String type, @Nullable final String treatment, @Nullable final String bestResponse) {
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
