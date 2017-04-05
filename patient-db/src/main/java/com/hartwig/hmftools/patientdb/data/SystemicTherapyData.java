package com.hartwig.hmftools.patientdb.data;

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

    public SystemicTherapyData(@Nullable final LocalDate startDate, @Nullable final LocalDate endDate,
            @Nullable final String type, @Nullable final String treatment, @Nullable final String bestResponse) {
        this.startDate = startDate;
        this.endDate = endDate;
        this.type = type;
        this.treatment = treatment;
        this.bestResponse = bestResponse;
    }

    @Nullable
    public LocalDate endDate() {
        return endDate;
    }

    @Nullable
    public LocalDate startDate() {
        return startDate;
    }

    @Nullable
    public String bestResponse() {
        return bestResponse;
    }

    @Nullable
    public String treatment() {
        return treatment;
    }

    @Nullable
    public String type() {
        return type;
    }
}
