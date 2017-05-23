package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientData;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.matchers.BiopsyMatcher;
import com.hartwig.hmftools.patientdb.matchers.TreatmentMatcher;
import com.hartwig.hmftools.patientdb.matchers.TreatmentResponseMatcher;

import org.apache.logging.log4j.ThreadContext;
import org.jetbrains.annotations.NotNull;

public class PatientReader {
    @NotNull
    private final CpctPatientReader cpctPatientReader;
    @NotNull
    private final SampleReader sampleReader;
    @NotNull
    private final BiopsyTreatmentReader biopsyTreatmentReader;

    public PatientReader(@NotNull final CpctEcrfModel model,
            @NotNull final Map<String, String> treatmentToTypeMappings, @NotNull final String limsCsv,
            @NotNull final String limsOldCsv, @NotNull final String umcuCsv) throws IOException, HartwigException {
        cpctPatientReader = new CpctPatientReader(model);
        sampleReader = new SampleReader(limsCsv, limsOldCsv, umcuCsv);
        biopsyTreatmentReader = new BiopsyTreatmentReader(treatmentToTypeMappings);
    }

    @NotNull
    public Patient read(@NotNull final EcrfPatient patient, @NotNull final List<String> tumorSamplesForPatient)
            throws IOException, HartwigException {
        ThreadContext.put("cpctHospitalCode", "HMF");
        final List<SampleData> sequencedBiopsies = sampleReader.read(tumorSamplesForPatient);
        ThreadContext.put("cpctHospitalCode", patient.patientId().substring(6, 8));
        final PatientData patientData = cpctPatientReader.read(patient);
        final List<BiopsyData> clinicalBiopsies = BiopsyReader.read(patient);
        final List<BiopsyTreatmentData> treatments = biopsyTreatmentReader.read(patient);
        final List<BiopsyTreatmentResponseData> treatmentResponses = BiopsyTreatmentResponseReader.read(patient);
        final List<BiopsyData> matchedBiopsies = BiopsyMatcher.matchBiopsiesToTumorSamples(patient.patientId(),
                sequencedBiopsies, clinicalBiopsies);
        final List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatmentsToBiopsies(
                patient.patientId(), clinicalBiopsies, treatments);
        final List<BiopsyTreatmentResponseData> matchedResponses = TreatmentResponseMatcher.matchTreatmentResponsesToTreatments(
                patient.patientId(), treatments, treatmentResponses);
        ThreadContext.put("cpctHospitalCode", "default");
        return new Patient(patientData, sequencedBiopsies, matchedBiopsies, matchedTreatments, matchedResponses);
    }
}
