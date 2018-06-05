package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Strings;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface TreatmentData extends Comparable<TreatmentData> {

    @VisibleForTesting
    String COMBI_THERAPY = "Combi therapy";

    @VisibleForTesting
    String COMBI_MECHANISM = "Combi mechanism";

    @Nullable
    String treatmentGiven();

    @Nullable
    String radiotherapyGiven();

    @NotNull
    List<DrugData> drugs();

    @NotNull
    FormStatus formStatus();

    @NotNull
    @Value.Derived
    default List<CuratedDrug> curatedDrugs() {
        return drugs().stream().flatMap(drug -> drug.filteredCuratedDrugs().stream()).collect(Collectors.toList());
    }

    @Nullable
    default String treatmentName() {
        List<CuratedDrug> drugs = curatedDrugs();
        Collections.sort(drugs);

        final String concatenatedTreatmentName = drugs.stream().map(CuratedDrug::name).collect(Collectors.joining("/"));
        return Strings.emptyToNull(concatenatedTreatmentName);
    }

    @Nullable
    default String concatenatedType() {
        List<CuratedDrug> drugs = curatedDrugs();
        Collections.sort(drugs);

        final String concatenatedTreatmentType = drugs.stream().map(CuratedDrug::type).collect(Collectors.joining("/"));
        return Strings.emptyToNull(concatenatedTreatmentType);
    }

    @Nullable
    default String consolidatedType() {
        final Set<String> types = curatedDrugs().stream().map(CuratedDrug::type).collect(Collectors.toSet());
        if (types.isEmpty()) {
            return null;
        } else if (types.size() == 1) {
            return types.iterator().next();
        } else {
            return COMBI_THERAPY;
        }
    }

    @Nullable
    default String concatenatedMechanism() {
        List<CuratedDrug> drugs = curatedDrugs();
        Collections.sort(drugs);

        final String concatenatedTreatmentMechanism = drugs.stream().map(CuratedDrug::mechanism).collect(Collectors.joining("/"));
        return Strings.emptyToNull(concatenatedTreatmentMechanism);
    }

    @Nullable
    default String consolidatedMechanism() {
        final Set<String> types = curatedDrugs().stream().map(CuratedDrug::mechanism).collect(Collectors.toSet());
        if (types.isEmpty()) {
            return null;
        } else if (types.size() == 1) {
            return types.iterator().next();
        } else {
            return COMBI_MECHANISM;
        }
    }

    @Nullable
    default LocalDate startDate() {
        LocalDate startDate = null;
        for (final DrugData drug : drugs()) {
            final LocalDate drugStartDate = drug.startDate();
            if (startDate == null || (drugStartDate != null && drugStartDate.isBefore(startDate))) {
                startDate = drugStartDate;
            }
        }
        return startDate;
    }

    @Nullable
    default LocalDate endDate() {
        if (drugs().isEmpty()) {
            return null;
        } else {
            LocalDate endDate = drugs().get(0).endDate();
            for (final DrugData drug : drugs()) {
                final LocalDate drugEndDate = drug.endDate();
                if (drugEndDate == null || (endDate != null && drugEndDate.isAfter(endDate))) {
                    endDate = drugEndDate;
                }
            }
            return endDate;
        }
    }

    @Override
    default int compareTo(@NotNull final TreatmentData other) {
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
