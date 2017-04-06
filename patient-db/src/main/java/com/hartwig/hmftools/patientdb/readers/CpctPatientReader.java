package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.consensus.ConsensusRule;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientInfo;
import com.hartwig.hmftools.patientdb.data.RadioTherapyData;
import com.hartwig.hmftools.patientdb.data.SomaticVariantData;
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
    @NotNull
    private final SomaticVariantReader somaticVariantReader;

    public CpctPatientReader(@NotNull final CpctEcrfModel model, @NotNull final ConsensusRule consensusRule)
            throws IOException, HartwigException {
        cpctPatientInfoReader = new CpctPatientInfoReader(model);
        cpctTumorDataReader = new CpctTumorDataReader();
        cpctSystemicTherapyReader = new CpctSystemicTherapyReader();
        cpctRadioTherapyReader = new CpctRadioTherapyReader();
        cpctTreatmentDataReader = new CpctTreatmentDataReader();
        somaticVariantReader = new SomaticVariantReader(consensusRule);
    }

    @NotNull
    public Patient read(@NotNull final EcrfPatient patient, @NotNull final String runDirectoryPath)
            throws IOException, HartwigException {
        ThreadContext.put("cpctHospitalCode", patient.patientId().substring(6, 8));
        final PatientInfo patientInfo = cpctPatientInfoReader.read(patient);
        final Optional<TumorData> tumorDataOpt = cpctTumorDataReader.read(patient);
        final Optional<List<SystemicTherapyData>> systemicTherapiesOpt = cpctSystemicTherapyReader.read(patient);
        final Optional<List<RadioTherapyData>> radioTherapiesOpt = cpctRadioTherapyReader.read(patient);
        final Optional<TreatmentData> treatmentDataOpt = cpctTreatmentDataReader.read(patient);
        final List<SomaticVariantData> somaticVariants = somaticVariantReader.read(runDirectoryPath);
        ThreadContext.put("cpctHospitalCode", "default");
        return new Patient(patientInfo, tumorDataOpt, systemicTherapiesOpt, radioTherapiesOpt, treatmentDataOpt,
                somaticVariants);
        //        return new Patient(patientInfo, tumorDataOpt, systemicTherapiesOpt, radioTherapiesOpt, treatmentDataOpt,
        //                Lists.newArrayList());
    }
}
