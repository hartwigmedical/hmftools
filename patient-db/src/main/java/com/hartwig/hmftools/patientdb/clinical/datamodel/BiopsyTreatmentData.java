package com.hartwig.hmftools.patientdb.clinical.datamodel;

import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyTreatmentData implements TreatmentData {

    private static final AtomicInteger ID_COUNTER = new AtomicInteger();

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

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    @NotNull
    public static BiopsyTreatmentData of(@Nullable Integer biopsyId, @Nullable String treatmentGiven, @Nullable String radiotherapyGiven,
            @NotNull List<DrugData> drugs, @NotNull FormStatus formStatus) {
        return ImmutableBiopsyTreatmentData.builder()
                .id(createId())
                .biopsyId(biopsyId)
                .treatmentGiven(treatmentGiven)
                .radiotherapyGiven(radiotherapyGiven)
                .drugs(drugs)
                .formStatus(formStatus)
                .build();
    }

    @Override
    public String toString() {
        return treatmentName() + " (" + startDate() + " - " + endDate() + ")";
    }
}
