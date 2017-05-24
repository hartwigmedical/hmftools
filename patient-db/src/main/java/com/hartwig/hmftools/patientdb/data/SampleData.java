package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.concurrent.atomic.AtomicInteger;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SampleData {

    @NotNull
    private final String sampleId;
    @NotNull
    private final LocalDate arrivalDate;
    @Nullable
    private final LocalDate samplingDate;

    public SampleData(@NotNull final String sampleId, @NotNull final LocalDate arrivalDate,
            @Nullable final LocalDate samplingDate) {
        this.sampleId = sampleId;
        this.arrivalDate = arrivalDate;
        this.samplingDate = samplingDate;
    }

    @NotNull
    public String sampleId() {
        return sampleId;
    }

    @NotNull
    public LocalDate arrivalDate() {
        return arrivalDate;
    }

    @Nullable
    public LocalDate samplingDate() {
        return samplingDate;
    }

    @NotNull
    public LocalDate date() {
        return samplingDate != null ? samplingDate : arrivalDate;
    }
}
