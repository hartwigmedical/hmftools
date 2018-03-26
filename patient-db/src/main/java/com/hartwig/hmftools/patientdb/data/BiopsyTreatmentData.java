package com.hartwig.hmftools.patientdb.data;

import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyTreatmentData implements TreatmentData {

    public abstract int id();

    @Nullable
    public abstract String treatmentGiven();

    @Nullable
    public abstract String radiotherapyGiven();

    @NotNull
    public abstract List<DrugData> drugs();

    @Nullable
    public abstract Integer biopsyId();

    @NotNull
    public abstract FormStatusState formStatus();

    public abstract boolean formLocked();

    private static final AtomicInteger ID_COUNTER = new AtomicInteger();

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    @NotNull
    public static BiopsyTreatmentData of(@Nullable final String treatmentGiven, @Nullable final String radiotherapyGiven,
            @NotNull final List<DrugData> drugs, @NotNull final FormStatusState formStatus, final boolean formLocked) {
        return ImmutableBiopsyTreatmentData.of(createId(), treatmentGiven, radiotherapyGiven, drugs, null, formStatus, formLocked);
    }

    @Override
    public String toString() {
        return treatmentName() + " (" + startDate() + " - " + endDate() + ")";
    }
}
