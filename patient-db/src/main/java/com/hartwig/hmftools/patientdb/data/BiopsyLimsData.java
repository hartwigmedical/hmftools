package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.concurrent.atomic.AtomicInteger;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BiopsyLimsData {

    private final int id;
    @NotNull
    private final String sampleId;
    @NotNull
    private final LocalDate arrivalDate;
    @Nullable
    private final LocalDate samplingDate;

    private static final AtomicInteger ID_COUNTER = new AtomicInteger();

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    public BiopsyLimsData(@NotNull final String sampleId, @NotNull final LocalDate arrivalDate,
            @Nullable final LocalDate samplingDate) {
        this.sampleId = sampleId;
        this.arrivalDate = arrivalDate;
        this.samplingDate = samplingDate;
        this.id = createId();
    }

    public int id() {
        return id;
    }

    @NotNull
    public LocalDate arrivalDate() {
        return arrivalDate;
    }

    @NotNull
    public String sampleId() {
        return sampleId;
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
