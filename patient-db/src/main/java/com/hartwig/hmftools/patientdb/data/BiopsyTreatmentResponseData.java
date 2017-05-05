package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import org.jetbrains.annotations.Nullable;

public class BiopsyTreatmentResponseData {

    @Nullable
    private final Integer treatmentId;

    @Nullable
    private final LocalDate assessmentDate;

    @Nullable
    private final LocalDate responseDate;

    @Nullable
    private final String response;

    @Nullable
    private final String measurementDone;

    public BiopsyTreatmentResponseData(@Nullable final LocalDate assessmentDate,
            @Nullable final LocalDate responseDate, @Nullable final String response,
            @Nullable final String measurementDone) {
        this(null, assessmentDate, responseDate, response, measurementDone);
    }

    public BiopsyTreatmentResponseData(@Nullable final Integer treatmentId, @Nullable final LocalDate assessmentDate,
            @Nullable final LocalDate responseDate, @Nullable final String response,
            @Nullable final String measurementDone) {
        this.treatmentId = treatmentId;
        this.assessmentDate = assessmentDate;
        this.responseDate = responseDate;
        this.response = response;
        this.measurementDone = measurementDone;
    }

    @Nullable
    public Integer treatmentId() {
        return treatmentId;
    }

    @Nullable
    public LocalDate date() {
        if (assessmentDate == null) {
            return responseDate;
        }
        return assessmentDate;
    }

    @Nullable
    public LocalDate assessmentDate() {
        return assessmentDate;
    }

    @Nullable
    public LocalDate responseDate() {
        return responseDate;
    }

    @Nullable
    public String response() {
        return response;
    }

    @Nullable
    public String measurementDone() {
        return measurementDone;
    }
}
