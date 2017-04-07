package com.hartwig.hmftools.patientdb.readers;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.TumorData;

import org.jetbrains.annotations.NotNull;

public class CpctTumorDataReader {
    private static final String FIELD_TUMORLOCATION = "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC";
    private static final String FIELD_ENTRYSTAGE = "BASELINE.CARCINOMA.CARCINOMA.ENTRYSTAGE";
    private static final String FIELD_BIOPSYLOCATION = "BIOPSY.BIOPS.BIOPSIES.BILESSITE";
    private static final String FIELD_BIOPSYLOCATIONOTHER = "BIOPSY.BIOPS.BIOPSIES.BIOTHLESSITE";

    @NotNull
    public Optional<TumorData> read(@NotNull final EcrfPatient patient) {
        final String tumorLocation = GenericReader.getField(patient, FIELD_TUMORLOCATION);
        final String tumorEntryStage = GenericReader.getField(patient, FIELD_ENTRYSTAGE);
        final List<String> biopsyLocationsField = GenericReader.getFieldValuesWithOthers(patient, FIELD_BIOPSYLOCATION,
                FIELD_BIOPSYLOCATIONOTHER);
        final List<String> biopsyLocations = biopsyLocationsField.stream().filter(
                location -> location != null && location.length() > 0).collect(Collectors.toList());
        if (Utils.anyNotNull(tumorLocation, tumorEntryStage) || biopsyLocations.size() > 0) {
            return Optional.of(new TumorData(tumorLocation, biopsyLocations, tumorEntryStage));
        } else {
            return Optional.empty();
        }
    }
}
