package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.concurrent.atomic.AtomicInteger;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyData implements Comparable<BiopsyData>  {

    public abstract int id();

    @Nullable
    public abstract String sampleId();

    @Nullable
    public abstract LocalDate date();

    @Nullable
    public abstract String biopsyTaken();

    @Nullable
    public abstract String biopsyEvaluable();

    @Nullable
    public abstract CuratedBiopsyType type();

    @Nullable
    public abstract String site();

    @Nullable
    public abstract String location();

    @NotNull
    public abstract FormStatus formStatus();

    private static final AtomicInteger ID_COUNTER = new AtomicInteger(0);

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    @NotNull
    public static BiopsyData of(@Nullable LocalDate date, @Nullable String biopsyTaken, @Nullable String biopsyEvaluable,
            @NotNull CuratedBiopsyType curatedBiopsyType, @Nullable String site, @Nullable String location,
            @NotNull FormStatus formStatus) {
        return ImmutableBiopsyData.of(createId(), null,
                date,
                biopsyTaken,
                biopsyEvaluable,
                curatedBiopsyType,
                site,
                location, formStatus);
    }

    @Nullable
    public String curatedType() {
        CuratedBiopsyType curatedType = type();
        return curatedType != null ? curatedType.type() : null;
    }

    public boolean isPotentiallyEvaluable() {
        String evaluable = biopsyEvaluable();
        return evaluable == null || !evaluable.equalsIgnoreCase("no");
    }

    @Override
    public String toString() {
        return site() + " - " + location() + " (" + date() + ")";
    }

    @Override
    public int compareTo(@NotNull final BiopsyData other) {
        LocalDate date1 = date();
        LocalDate date2 = other.date();
        if (date1 == null && date2 == null) {
            return 0;
        } else if (date1 == null) {
            return 1;
        } else if (date2 == null) {
            return -1;
        } else {
            return date1.compareTo(date2);
        }
    }
}
