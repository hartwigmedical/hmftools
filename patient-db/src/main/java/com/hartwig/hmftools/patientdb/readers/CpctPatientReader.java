package com.hartwig.hmftools.patientdb.readers;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientInfo;
import com.hartwig.hmftools.patientdb.data.RadioTherapyData;
import com.hartwig.hmftools.patientdb.data.SystemicTherapyData;
import com.hartwig.hmftools.patientdb.data.TreatmentData;
import com.hartwig.hmftools.patientdb.data.TumorData;

import org.apache.logging.log4j.ThreadContext;
import org.jetbrains.annotations.NotNull;

public class CpctPatientReader {
    @NotNull
    private final CpctPatientInfoReader cpctPatientInfoReader;
    @NotNull
    private final CpctTumorDataReader cpctTumorDataReader;
    @NotNull
    private final CpctSystemicTherapyReader cpctSystemicTherapyReader;
    @NotNull
    private final CpctRadioTherapyReader cpctRadioTherapyReader;
    @NotNull
    private final CpctTreatmentDataReader cpctTreatmentDataReader;

    public CpctPatientReader(@NotNull final CpctEcrfModel model) {
        cpctPatientInfoReader = new CpctPatientInfoReader(model);
        cpctTumorDataReader = new CpctTumorDataReader();
        cpctSystemicTherapyReader = new CpctSystemicTherapyReader();
        cpctRadioTherapyReader = new CpctRadioTherapyReader();
        cpctTreatmentDataReader = new CpctTreatmentDataReader();
    }

    @NotNull
    public Patient read(@NotNull final EcrfPatient patient) {
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
