package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.consensus.ConsensusRule;
import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.BiopsyLimsData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientInfo;
import com.hartwig.hmftools.patientdb.matchers.BiopsyMatcher;
import com.hartwig.hmftools.patientdb.matchers.TreatmentMatcher;
import com.hartwig.hmftools.patientdb.matchers.TreatmentResponseMatcher;

import org.apache.logging.log4j.ThreadContext;
import org.jetbrains.annotations.NotNull;

public class CpctPatientReader {
    @NotNull
    private final CpctPatientInfoReader cpctPatientInfoReader;
    @NotNull
    private final BiopsyLimsDataReader biopsyLimsDataReader;
    @NotNull
    private final BiopsyTreatmentReader biopsyTreatmentReader;
    @NotNull
    private final SomaticVariantReader somaticVariantReader;

    public CpctPatientReader(@NotNull final CpctEcrfModel model, @NotNull final ConsensusRule consensusRule,
            @NotNull final Map<String, String> treatmentMappings, @NotNull final String limsCsv,
            @NotNull final String limsOldCsv) throws IOException, HartwigException {
        cpctPatientInfoReader = new CpctPatientInfoReader(model);
        biopsyLimsDataReader = new BiopsyLimsDataReader(limsCsv, limsOldCsv);
        biopsyTreatmentReader = new BiopsyTreatmentReader(treatmentMappings);
        somaticVariantReader = new SomaticVariantReader(consensusRule);
    }

    @NotNull
    public Patient read(@NotNull final EcrfPatient patient, @NotNull final List<String> sampleIdsForPatient)
            throws IOException, HartwigException {
        ThreadContext.put("cpctHospitalCode", patient.patientId().substring(6, 8));
        final PatientInfo patientInfo = cpctPatientInfoReader.read(patient);
        final List<BiopsyClinicalData> clinicalBiopsies = BiopsyClinicalDataReader.read(patient);
        final List<BiopsyLimsData> sequencedBiopsies = biopsyLimsDataReader.read(sampleIdsForPatient);
        final List<BiopsyTreatmentData> treatments = biopsyTreatmentReader.read(patient);
        final List<BiopsyTreatmentResponseData> treatmentResponses = BiopsyTreatmentResponseReader.read(patient);
        final List<BiopsyTreatmentResponseData> matchedResponses = TreatmentResponseMatcher.matchTreatmentResponses(
                patient.patientId(), treatments, treatmentResponses);
        final List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatments(patient.patientId(),
                clinicalBiopsies, treatments);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies(patient.patientId(),
                sequencedBiopsies, clinicalBiopsies);
        //MIVO: skip reading somatic variants for now
        //        final List<SomaticVariantData> somaticVariants = somaticVariantReader.read(runDirectoryPath);
        ThreadContext.put("cpctHospitalCode", "default");
        //        return new Patient(patientInfo, tumorDataOpt, systemicTherapies, radioTherapies, treatmentDataOpt,
        //                somaticVariants);
        return new Patient(patientInfo, sequencedBiopsies, matchedTreatments, matchedResponses, matchedBiopsies,
                Lists.newArrayList());
    }

}
