package com.hartwig.hmftools.patientdb.readers;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.data.TumorData;

import org.jetbrains.annotations.NotNull;

class CpctTumorDataReader {
    private static final String FIELD_TUMORLOCATION = "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC";
    private static final String FIELD_ENTRYSTAGE = "BASELINE.CARCINOMA.CARCINOMA.ENTRYSTAGE";
    private static final String FIELD_BIOPSYLOCATION = "BIOPSY.BIOPS.BIOPSIES.BILESSITE";
    private static final String FIELD_BIOPSYLOCATIONOTHER = "BIOPSY.BIOPS.BIOPSIES.BIOTHLESSITE";

    @NotNull
    Optional<TumorData> read(@NotNull final EcrfPatient patient) {
        final String tumorLocation = GenericReader.getField(patient, FIELD_TUMORLOCATION);
        final String tumorEntryStage = GenericReader.getField(patient, FIELD_ENTRYSTAGE);
        final List<String> biopsyLocations = GenericReader.getFieldValuesWithOthers(patient, FIELD_BIOPSYLOCATION,
                FIELD_BIOPSYLOCATIONOTHER);
        if ((tumorLocation == null || tumorLocation.replaceAll("\\s", "").length() == 0) && (tumorEntryStage == null
                || tumorEntryStage.replaceAll("\\s", "").length() == 0) && (biopsyLocations == null
                || biopsyLocations.size() == 0)) {
            return Optional.empty();
        } else {
            return Optional.of(new TumorData(tumorLocation, biopsyLocations, tumorEntryStage));
        }
    }
}
