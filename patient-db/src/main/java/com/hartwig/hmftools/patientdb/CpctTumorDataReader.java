package com.hartwig.hmftools.patientdb;

import java.util.List;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.jetbrains.annotations.NotNull;

public class CpctTumorDataReader {
    @NotNull
    public TumorData read(@NotNull EcrfPatient patient) {
        final String tumorLocation = GenericReader.getField(patient, "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC");
        final String tumorEntryStage = GenericReader.getField(patient, "BASELINE.CARCINOMA.CARCINOMA.ENTRYSTAGE");
        final List<String> biopsyLocations = GenericReader.getFieldValues(patient, "BIOPSY.BIOPS.BIOPSIES.BILESSITE");
        return new TumorData(tumorLocation, biopsyLocations, tumorEntryStage);
    }
}
