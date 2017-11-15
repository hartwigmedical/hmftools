package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Strings;
import com.hartwig.hmftools.patientdb.Utils;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyTreatmentData {

    @VisibleForTesting
    static final String COMBI_THERAPY = "Combi therapy";

    public abstract int id();

    @Nullable
    public abstract String treatmentGiven();

    @Nullable
    public abstract LocalDate startDate();

    @Nullable
    public abstract LocalDate endDate();

    @NotNull
    public abstract List<BiopsyTreatmentDrugData> drugs();

    @Nullable
    public abstract Integer biopsyId();

    @NotNull
    public abstract String formStatus();

    @NotNull
    public abstract String formLocked();

    @Value.Derived
    List<CuratedTreatment> curatedDrugs() {
        return drugs().stream().flatMap(drug -> drug.filteredCuratedTreatments().stream()).collect(Collectors.toList());
    }

    private static final AtomicInteger ID_COUNTER = new AtomicInteger();

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    @NotNull
    public static BiopsyTreatmentData of(@Nullable final String treatmentGiven, @Nullable final LocalDate startDate,
            @Nullable final LocalDate endDate, @NotNull final List<BiopsyTreatmentDrugData> drugs, @NotNull final String formStatus,
            @NotNull final String formLocked) {

        return ImmutableBiopsyTreatmentData.of(createId(), treatmentGiven, startDate, endDate, drugs, null, formStatus, formLocked);
    }

    @Nullable
    public String treatmentName() {
        final String distinctSortedDrugs = curatedDrugs().stream()
                .map(treatment -> Utils.capitalize(treatment.name()))
                .sorted()
                .distinct()
                .collect(Collectors.joining("/"));
        return Strings.emptyToNull(distinctSortedDrugs);
    }

    @Nullable
    public String type() {
        final Set<String> types = curatedDrugs().stream().map(CuratedTreatment::type).collect(Collectors.toSet());
        if (types.isEmpty()) {
            return null;
        } else if (types.size() == 1) {
            return types.iterator().next();
        } else {
            return COMBI_THERAPY;
        }
    }

    @Override
    public String toString() {
        return treatmentName() + "(" + startDate() + " - " + endDate() + ")";
    }
}
