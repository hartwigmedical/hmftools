package com.hartwig.hmftools.patientdb;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.jetbrains.annotations.NotNull;

class CpctPatientReader {
    private final CpctPatientInfoReader cpctPatientInfoReader;
    private final CpctTumorDataReader cpctTumorDataReader;
    private final CpctSystemicTherapyReader cpctSystemicTherapyReader;
    private final CpctRadioTherapyReader cpctRadioTherapyReader;
    private final CpctTreatmentDataReader cpctTreatmentDataReader;

    CpctPatientReader(@NotNull CpctEcrfModel model) {
        cpctPatientInfoReader = new CpctPatientInfoReader(model);
        cpctTumorDataReader = new CpctTumorDataReader();
        cpctSystemicTherapyReader = new CpctSystemicTherapyReader();
        cpctRadioTherapyReader = new CpctRadioTherapyReader();
        cpctTreatmentDataReader = new CpctTreatmentDataReader();
    }

    @NotNull
    Patient read(@NotNull EcrfPatient patient) {
        final PatientInfo patientInfo = cpctPatientInfoReader.read(patient);
        final Optional<TumorData> tumorDataOpt = cpctTumorDataReader.read(patient);
        final Optional<List<SystemicTherapyData>> systemicTherapiesOpt = cpctSystemicTherapyReader.read(patient);
        final Optional<List<RadioTherapyData>> radioTherapiesOpt = cpctRadioTherapyReader.read(patient);
        final Optional<TreatmentData> treatmentDataOpt = cpctTreatmentDataReader.read(patient);
        return new Patient(patientInfo, tumorDataOpt, systemicTherapiesOpt, radioTherapiesOpt, treatmentDataOpt);
    }
}
