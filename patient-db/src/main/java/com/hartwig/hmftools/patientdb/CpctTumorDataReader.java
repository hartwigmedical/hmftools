package com.hartwig.hmftools.patientdb;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.jetbrains.annotations.NotNull;

class CpctTumorDataReader {
    @NotNull
    Optional<TumorData> read(@NotNull EcrfPatient patient) {
        final String tumorLocation = GenericReader.getField(patient, "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC");
        final String tumorEntryStage = GenericReader.getField(patient, "BASELINE.CARCINOMA.CARCINOMA.ENTRYSTAGE");
        final List<String> biopsyLocations = GenericReader.getFieldValues(patient, "BIOPSY.BIOPS.BIOPSIES.BILESSITE");
        if ((tumorLocation == null || tumorLocation.replaceAll("\\s", "").length() == 0) && (tumorEntryStage == null
                || tumorEntryStage.replaceAll("\\s", "").length() == 0) && (biopsyLocations == null
                || biopsyLocations.size() == 0)) {
            return Optional.empty();
        } else {
            return Optional.of(new TumorData(tumorLocation, biopsyLocations, tumorEntryStage));
        }
    }
}
