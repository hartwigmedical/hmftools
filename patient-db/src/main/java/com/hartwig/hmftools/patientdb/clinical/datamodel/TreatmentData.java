package com.hartwig.hmftools.patientdb.clinical.datamodel;

import java.time.LocalDate;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Strings;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface TreatmentData extends Comparable<TreatmentData> {

    String SEPARATOR = "/";

    @VisibleForTesting
    String COMBI_TYPE = "Multiple therapy";

    @VisibleForTesting
    String COMBI_MECHANISM = "Multiple mechanism";

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
    @Value.Derived
    default String treatmentName() {
        List<CuratedDrug> drugs = curatedDrugs();

        String concatenatedTreatmentName = drugs.stream().map(CuratedDrug::name).collect(Collectors.joining(SEPARATOR));
        return Strings.emptyToNull(concatenatedTreatmentName);
    }

    @Nullable
    @Value.Derived
    default String concatenatedType() {
        List<CuratedDrug> drugs = curatedDrugs();

        String value = drugs.stream().map(CuratedDrug::type).collect(Collectors.joining(SEPARATOR));
        return Strings.emptyToNull(value);
    }

    @Nullable
    @Value.Derived
    default String concatenatedMechanism() {
        List<CuratedDrug> drugs = curatedDrugs();

        String value = drugs.stream().map(CuratedDrug::mechanism).collect(Collectors.joining(SEPARATOR));
        return Strings.emptyToNull(value);
    }

    @Nullable
    @Value.Derived
    default String consolidatedType() {
        return consolidate(curatedDrugs().stream().map(CuratedDrug::type).collect(Collectors.toSet()), COMBI_TYPE);
    }

    @Nullable
    @Value.Derived
    default String consolidatedMechanism() {
        return consolidate(curatedDrugs().stream().map(CuratedDrug::mechanism).collect(Collectors.toSet()), COMBI_MECHANISM);
    }

    @Nullable
    static String consolidate(@NotNull Set<String> values, @NotNull String combiValue) {
        if (values.isEmpty()) {
            return null;
        } else if (values.size() == 1) {
            return values.iterator().next();
        } else {
            return combiValue;
        }
    }

    @Nullable
    @Value.Derived
    default LocalDate startDate() {
        LocalDate startDate = null;
        for (DrugData drug : drugs()) {
            LocalDate drugStartDate = drug.startDate();
            if (startDate == null || (drugStartDate != null && drugStartDate.isBefore(startDate))) {
                startDate = drugStartDate;
            }
        }
        return startDate;
    }

    @Nullable
    @Value.Derived
    default LocalDate endDate() {
        if (drugs().isEmpty()) {
            return null;
        } else {
            LocalDate endDate = drugs().get(0).endDate();
            for (DrugData drug : drugs()) {
                LocalDate drugEndDate = drug.endDate();
                if (drugEndDate == null || (endDate != null && drugEndDate.isAfter(endDate))) {
                    endDate = drugEndDate;
                }
            }
            return endDate;
        }
    }

    @Override
    default int compareTo(@NotNull TreatmentData other) {
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
