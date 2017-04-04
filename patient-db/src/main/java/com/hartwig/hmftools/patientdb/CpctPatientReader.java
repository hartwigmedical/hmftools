package com.hartwig.hmftools.patientdb;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.ThreadContext;
import org.jetbrains.annotations.NotNull;

public class CpctPatientReader {
    private final CpctPatientInfoReader cpctPatientInfoReader;
    private final CpctTumorDataReader cpctTumorDataReader;
    private final CpctSystemicTherapyReader cpctSystemicTherapyReader;
    private final CpctRadioTherapyReader cpctRadioTherapyReader;
    private final CpctTreatmentDataReader cpctTreatmentDataReader;

    public CpctPatientReader(@NotNull CpctEcrfModel model) {
        cpctPatientInfoReader = new CpctPatientInfoReader(model);
        cpctTumorDataReader = new CpctTumorDataReader();
        cpctSystemicTherapyReader = new CpctSystemicTherapyReader();
        cpctRadioTherapyReader = new CpctRadioTherapyReader();
        cpctTreatmentDataReader = new CpctTreatmentDataReader();
    }

    @NotNull
    public Patient read(@NotNull EcrfPatient patient) {
        ThreadContext.put("cpctHospitalCode", patient.patientId().substring(6, 8));
        final PatientInfo patientInfo = cpctPatientInfoReader.read(patient);
        final Optional<TumorData> tumorDataOpt = cpctTumorDataReader.read(patient);
        final Optional<List<SystemicTherapyData>> systemicTherapiesOpt = cpctSystemicTherapyReader.read(patient);
        final Optional<List<RadioTherapyData>> radioTherapiesOpt = cpctRadioTherapyReader.read(patient);
        final Optional<TreatmentData> treatmentDataOpt = cpctTreatmentDataReader.read(patient);
        ThreadContext.put("cpctHospitalCode", "default");
        return new Patient(patientInfo, tumorDataOpt, systemicTherapiesOpt, radioTherapiesOpt, treatmentDataOpt);
    }
}
