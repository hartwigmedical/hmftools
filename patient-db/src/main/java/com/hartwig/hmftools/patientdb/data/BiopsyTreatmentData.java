package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Strings;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.Utils;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyTreatmentData implements Comparable<BiopsyTreatmentData> {

    @VisibleForTesting
    static final String COMBI_THERAPY = "Combi therapy";

    public abstract int id();

    @Nullable
    public abstract String treatmentGiven();

    @NotNull
    public abstract List<BiopsyTreatmentDrugData> drugs();

    @Nullable
    public abstract Integer biopsyId();

    @NotNull
    public abstract FormStatusState formStatus();

    public abstract boolean formLocked();

    @Value.Derived
    List<CuratedTreatment> curatedDrugs() {
        return drugs().stream().flatMap(drug -> drug.filteredCuratedTreatments().stream()).collect(Collectors.toList());
    }

    private static final AtomicInteger ID_COUNTER = new AtomicInteger();

    private static int createId() {
        return ID_COUNTER.getAndIncrement();
    }

    @NotNull
    public static BiopsyTreatmentData of(@Nullable final String treatmentGiven, @NotNull final List<BiopsyTreatmentDrugData> drugs,
            @NotNull final FormStatusState formStatus, final boolean formLocked) {
        return ImmutableBiopsyTreatmentData.of(createId(), treatmentGiven, drugs, null, formStatus, formLocked);
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

    @Nullable
    public LocalDate startDate() {
        LocalDate startDate = null;
        for (final BiopsyTreatmentDrugData drug : drugs()) {
            final LocalDate drugStartDate = drug.startDate();
            if (startDate == null || (drugStartDate != null && drugStartDate.isBefore(startDate))) {
                startDate = drugStartDate;
            }
        }
        return startDate;
    }

    @Nullable
    public LocalDate endDate() {
        if (drugs().isEmpty()) {
            return null;
        } else {
            LocalDate endDate = drugs().get(0).endDate();
            for (final BiopsyTreatmentDrugData drug : drugs()) {
                final LocalDate drugEndDate = drug.endDate();
                if (drugEndDate == null || (endDate != null && drugEndDate.isAfter(endDate))) {
                    endDate = drugEndDate;
                }
            }
            return endDate;
        }
    }

    @Override
    public String toString() {
        return treatmentName() + "(" + startDate() + " - " + endDate() + ")";
    }

    @Override
    public int compareTo(@NotNull final BiopsyTreatmentData other) {
        LocalDate startDate1 = startDate();
        LocalDate startDate2 = other.startDate();
        if (startDate1 == null && startDate2 == null) {
            return 0;
        } else if (startDate1 == null) {
            return 1;
        } else if (startDate2 == null) {
            return -1;
        } else {
            return startDate1.compareTo(startDate2);
        }
    }
}
