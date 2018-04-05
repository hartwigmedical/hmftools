package com.hartwig.hmftools.patientdb.data;

import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyTreatmentData implements TreatmentData {

    public abstract int id();

    @Nullable
    public abstract Integer biopsyId();

    @Nullable
    @Override
    public abstract String treatmentGiven();

    @Nullable
    @Override
    public abstract String radiotherapyGiven();

    @NotNull
    @Override
    public abstract List<DrugData> drugs();

    @NotNull
    @Override
    public abstract FormStatus formStatus();

    private static final AtomicInteger ID_COUNTER = new AtomicInteger();

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    @NotNull
    public static BiopsyTreatmentData of(@Nullable final String treatmentGiven, @Nullable final String radiotherapyGiven,
            @NotNull final List<DrugData> drugs, @NotNull final FormStatus formStatus) {
        return ImmutableBiopsyTreatmentData.of(createId(), null, treatmentGiven, radiotherapyGiven, drugs, formStatus);
    }

    @Override
    public String toString() {
        return treatmentName() + " (" + startDate() + " - " + endDate() + ")";
    }
}
