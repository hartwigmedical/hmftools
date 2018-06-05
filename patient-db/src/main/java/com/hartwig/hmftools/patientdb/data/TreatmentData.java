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
    String seperator = "/";

    @NotNull
    @Value.Derived
    default List<CuratedDrug> curatedDrugs() {
        return drugs().stream().flatMap(drug -> drug.filteredCuratedDrugs().stream()).collect(Collectors.toList());
    }

    @Nullable
    default String treatmentName() {
        List<CuratedDrug> drugs = curatedDrugs();
        Collections.sort(drugs);

        final String concatenatedTreatmentName = drugs.stream().map(CuratedDrug::name).collect(Collectors.joining(seperator));
        return Strings.emptyToNull(concatenatedTreatmentName);
    }

    @Nullable
    default String concatenatedTypeOrMechanism(@NotNull String ValueForTreatment) {
        List<CuratedDrug> drugs = curatedDrugs();
        Collections.sort(drugs);

        return Strings.emptyToNull(chcekValueConcatenated(ValueForTreatment, drugs));
    }

    @Nullable
    default String chcekValueConcatenated(@NotNull String TreatmentValue, @NotNull List<CuratedDrug> drugs) {
        return TreatmentValue.equals("type")
                ? drugs.stream().map(CuratedDrug::type).collect(Collectors.joining(seperator))
                : drugs.stream().map(CuratedDrug::mechanism).collect(Collectors.joining(seperator));
    }

    @Nullable
    default Set<String> checkValueConsolidated(@NotNull String TreatmentValue) {
        return TreatmentValue.equals("type")
                ? curatedDrugs().stream().map(CuratedDrug::type).collect(Collectors.toSet())
                : curatedDrugs().stream().map(CuratedDrug::mechanism).collect(Collectors.toSet());
    }

    @NotNull
    default String checkCombiValueTreatment(@NotNull String valueCombi) {
        return valueCombi.equals("type") ? COMBI_THERAPY : COMBI_MECHANISM;
    }

    @Nullable
    default String consolidatedTypeOrMechanism(@NotNull String ValueForTreatment) {
        if (checkValueConsolidated(ValueForTreatment).isEmpty()) {
            return null;
        } else if (checkValueConsolidated(ValueForTreatment).size() == 1) {
            return checkValueConsolidated(ValueForTreatment).iterator().next();
        } else {
            return checkCombiValueTreatment(ValueForTreatment);
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
