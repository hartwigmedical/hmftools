package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.consensus.ConsensusRule;
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
    @NotNull
    private final SomaticVariantReader somaticVariantReader;

    public CpctPatientReader(@NotNull final CpctEcrfModel model, @NotNull final ConsensusRule consensusRule,
            @NotNull final String treatmentMappingCsv) throws IOException, HartwigException {
        cpctPatientInfoReader = new CpctPatientInfoReader(model);
        cpctTumorDataReader = new CpctTumorDataReader();
        cpctSystemicTherapyReader = new CpctSystemicTherapyReader();
        cpctRadioTherapyReader = new CpctRadioTherapyReader();
        cpctTreatmentDataReader = new CpctTreatmentDataReader(treatmentMappingCsv);
        somaticVariantReader = new SomaticVariantReader(consensusRule);
    }

    @NotNull
    public Patient read(@NotNull final EcrfPatient patient, @NotNull final String runDirectoryPath)
            throws IOException, HartwigException {
        ThreadContext.put("cpctHospitalCode", patient.patientId().substring(6, 8));
        final PatientInfo patientInfo = cpctPatientInfoReader.read(patient);
        final Optional<TumorData> tumorDataOpt = cpctTumorDataReader.read(patient);
        final List<SystemicTherapyData> systemicTherapies = cpctSystemicTherapyReader.read(patient);
        final List<RadioTherapyData> radioTherapies = cpctRadioTherapyReader.read(patient);
        final Optional<TreatmentData> treatmentDataOpt = cpctTreatmentDataReader.read(patient);
        //MIVO: skip reading somatic variants for now
        //        final List<SomaticVariantData> somaticVariants = somaticVariantReader.read(runDirectoryPath);
        ThreadContext.put("cpctHospitalCode", "default");
        //        return new Patient(patientInfo, tumorDataOpt, systemicTherapies, radioTherapies, treatmentDataOpt,
        //                somaticVariants);
        return new Patient(patientInfo, tumorDataOpt, systemicTherapies, radioTherapies, treatmentDataOpt,
                Lists.newArrayList());
    }
}
