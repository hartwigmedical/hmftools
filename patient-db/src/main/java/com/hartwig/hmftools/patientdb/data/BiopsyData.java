package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.concurrent.atomic.AtomicInteger;

import org.jetbrains.annotations.Nullable;

public class BiopsyData {

    private final int id;
    @Nullable
    private final LocalDate date;
    @Nullable
    private final String location;
    @Nullable
    private final String sampleId;

    private static final AtomicInteger ID_COUNTER = new AtomicInteger(0);

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    public BiopsyData(@Nullable final LocalDate date, @Nullable final String location) {
        this(createId(), date, location, null);
    }

    public BiopsyData(final int id, @Nullable final LocalDate date, @Nullable final String location,
            @Nullable final String sampleId) {
        this.id = id;
        this.date = date;
        this.location = location;
        this.sampleId = sampleId;
    }

    public int id() {
        return id;
    }

    @Nullable
    public LocalDate date() {
        return date;
    }

    @Nullable
    public String location() {
        return location;
    }

    @Nullable
    public String sampleId() {
        return sampleId;
    }
}
