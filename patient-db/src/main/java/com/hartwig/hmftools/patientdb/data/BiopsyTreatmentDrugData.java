package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyTreatmentDrugData {
    @Nullable
    public abstract String name();

    @Nullable
    public abstract String type();

    @Nullable
    public abstract LocalDate startDate();

    @Nullable
    public abstract LocalDate endDate();
}
