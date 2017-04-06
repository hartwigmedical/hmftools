package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import org.jetbrains.annotations.Nullable;

public class TreatmentData {
    //    private final String type;
    @Nullable
    private final String treatmentName;
    @Nullable
    private final LocalDate startDate;
    @Nullable
    private final LocalDate endDate;
    @Nullable
    private final String earlyResponse;
    @Nullable
    private final LocalDate radiotherapyStartDate;
    @Nullable
    private final LocalDate radiotherapyEndDate;

    public TreatmentData(@Nullable final LocalDate startDate, @Nullable final LocalDate endDate,
            @Nullable final String treatmentName, @Nullable final String earlyResponse,
            @Nullable final LocalDate radiotherapyStartDate, @Nullable final LocalDate radiotherapyEndDate) {
        this.startDate = startDate;
        this.endDate = endDate;
        this.treatmentName = treatmentName;
        this.earlyResponse = earlyResponse;
        this.radiotherapyStartDate = radiotherapyStartDate;
        this.radiotherapyEndDate = radiotherapyEndDate;
    }

    @Nullable
    public String treatmentName() {
        return treatmentName;
    }

    @Nullable
    public LocalDate startDate() {
        return startDate;
    }

    @Nullable
    public LocalDate endDate() {
        return endDate;
    }

    @Nullable
    public String earlyResponse() {
        return earlyResponse;
    }

    @Nullable
    public LocalDate radiotherapyStartDate() {
        return radiotherapyStartDate;
    }

    @Nullable
    public LocalDate radiotherapyEndDate() {
        return radiotherapyEndDate;
    }
}
