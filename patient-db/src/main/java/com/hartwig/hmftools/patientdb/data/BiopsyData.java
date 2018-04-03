package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.concurrent.atomic.AtomicInteger;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

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
    public abstract String biopsyTaken();

    @Nullable
    public abstract String biopsyEvaluable();

    @Nullable
    public abstract CuratedBiopsyType type();

    @Nullable
    public String curatedType() {
        CuratedBiopsyType curatedType = type();
        return curatedType != null ? curatedType.type() : null;
    }

    @Nullable
    public abstract String site();

    @Nullable
    public abstract String location();

    @Nullable
    public abstract String sampleId();

    @NotNull
    public abstract FormStatusState formStatus();

    public abstract boolean formLocked();

    private static final AtomicInteger ID_COUNTER = new AtomicInteger(0);

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    @NotNull
    public static BiopsyData of(@Nullable LocalDate date, @Nullable String biopsyTaken, @Nullable String biopsyEvaluable,
            @NotNull CuratedBiopsyType curatedBiopsyType, @Nullable String site, @Nullable String location,
            @NotNull FormStatusState formStatus, boolean formLocked) {
        return ImmutableBiopsyData.of(createId(),
                date,
                biopsyTaken,
                biopsyEvaluable,
                curatedBiopsyType,
                site,
                location,
                null,
                formStatus,
                formLocked);
    }

    public boolean isPotentiallyEvaluable() {
        String evaluable = biopsyEvaluable();
        return evaluable == null || !evaluable.equalsIgnoreCase("no");
    }

    @Override
    public String toString() {
        return site() + " - " + location() + " (" + date() + ")";
    }
}
