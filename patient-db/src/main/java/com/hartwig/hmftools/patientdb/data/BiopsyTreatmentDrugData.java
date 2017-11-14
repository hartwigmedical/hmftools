package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.List;
import java.util.stream.Collectors;

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
    public abstract LocalDate startDate();

    @Nullable
    public abstract LocalDate endDate();

    @NotNull
    public abstract List<CuratedTreatment> curatedTreatments();

    @NotNull
    @Value.Derived
    public List<CuratedTreatment> filteredCuratedTreatments() {
        return curatedTreatments().stream()
                .filter(curatedTreatment -> !curatedTreatment.type().toLowerCase().equals("remove"))
                .distinct()
                .collect(Collectors.toList());
    }
}
