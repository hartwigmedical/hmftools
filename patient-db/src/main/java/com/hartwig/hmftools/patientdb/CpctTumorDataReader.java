package com.hartwig.hmftools.patientdb;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.jetbrains.annotations.NotNull;

class CpctTumorDataReader {
    private static final String FIELD_TUMORLOCATION = "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC";
    private static final String FIELD_ENTRYSTAGE = "BASELINE.CARCINOMA.CARCINOMA.ENTRYSTAGE";
    private static final String FIELD_BIOPSYLOCATION = "BIOPSY.BIOPS.BIOPSIES.BILESSITE";

    @NotNull
    Optional<TumorData> read(@NotNull EcrfPatient patient) {
        final String tumorLocation = GenericReader.getField(patient, FIELD_TUMORLOCATION);
        final String tumorEntryStage = GenericReader.getField(patient, FIELD_ENTRYSTAGE);
        final List<String> biopsyLocations = GenericReader.getFieldValues(patient, FIELD_BIOPSYLOCATION);
        if ((tumorLocation == null || tumorLocation.replaceAll("\\s", "").length() == 0) && (tumorEntryStage == null
                || tumorEntryStage.replaceAll("\\s", "").length() == 0) && (biopsyLocations == null
                || biopsyLocations.size() == 0)) {
            return Optional.empty();
        } else {
            return Optional.of(new TumorData(tumorLocation, biopsyLocations, tumorEntryStage));
        }
    }
}
