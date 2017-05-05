package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.concurrent.atomic.AtomicInteger;

import org.jetbrains.annotations.NotNull;

public class BiopsyLimsData {

    private final int id;
    @NotNull
    private final String sampleId;
    @NotNull
    private final LocalDate arrivalDate;

    private static final AtomicInteger idCounter = new AtomicInteger();

    private static int createId() {
        return idCounter.getAndIncrement();
    }

    public BiopsyLimsData(@NotNull final String sampleId, @NotNull final LocalDate arrivalDate) {
        this.sampleId = sampleId;
        this.arrivalDate = arrivalDate;
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
}
