package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.concurrent.atomic.AtomicInteger;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyData {

    public abstract int id();

    @Nullable
    public abstract LocalDate date();

    @Nullable
    public abstract String location();

    @Nullable
    public abstract String sampleId();

    @NotNull
    public abstract String formStatus();

    @NotNull
    public abstract String formLocked();

    private static final AtomicInteger ID_COUNTER = new AtomicInteger(0);

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    @NotNull
    public static BiopsyData of(@Nullable final LocalDate date, @Nullable final String location, @NotNull final String formStatus,
            @NotNull final String formLocked) {
        return ImmutableBiopsyData.of(createId(), date, location, null, formStatus, formLocked);
    }

    @Override
    public String toString() {
        return location() + "(" + date() + ")";
    }
}
