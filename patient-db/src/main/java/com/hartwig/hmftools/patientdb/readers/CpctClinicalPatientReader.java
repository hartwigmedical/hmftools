package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.exception.HartwigException;
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

public class CpctClinicalPatientReader {
    @NotNull
    private final CpctPatientInfoReader cpctPatientInfoReader;
    @NotNull
    private final BiopsyLimsDataReader biopsyLimsDataReader;
    @NotNull
    private final BiopsyTreatmentReader biopsyTreatmentReader;

    public CpctClinicalPatientReader(@NotNull final CpctEcrfModel model,
            @NotNull final Map<String, String> treatmentToTypeMappings, @NotNull final String limsCsv,
            @NotNull final String limsOldCsv, @NotNull final String umcuCsv) throws IOException, HartwigException {
        cpctPatientInfoReader = new CpctPatientInfoReader(model);
        biopsyLimsDataReader = new BiopsyLimsDataReader(limsCsv, limsOldCsv, umcuCsv);
        biopsyTreatmentReader = new BiopsyTreatmentReader(treatmentToTypeMappings);
    }

    @NotNull
    public Patient read(@NotNull final EcrfPatient patient, @NotNull final List<String> sampleIdsForPatient)
            throws IOException, HartwigException {
        ThreadContext.put("cpctHospitalCode", "HMF");
        final List<BiopsyLimsData> sequencedBiopsies = biopsyLimsDataReader.read(sampleIdsForPatient);
        ThreadContext.put("cpctHospitalCode", patient.patientId().substring(6, 8));
        final PatientInfo patientInfo = cpctPatientInfoReader.read(patient);
        final List<BiopsyClinicalData> clinicalBiopsies = BiopsyClinicalDataReader.read(patient);
        final List<BiopsyTreatmentData> treatments = biopsyTreatmentReader.read(patient);
        final List<BiopsyTreatmentResponseData> treatmentResponses = BiopsyTreatmentResponseReader.read(patient);
        final List<BiopsyTreatmentResponseData> matchedResponses = TreatmentResponseMatcher.matchTreatmentResponses(
                patient.patientId(), treatments, treatmentResponses);
        final List<BiopsyTreatmentData> matchedTreatments = TreatmentMatcher.matchTreatments(patient.patientId(),
                clinicalBiopsies, treatments);
        final List<BiopsyClinicalData> matchedBiopsies = BiopsyMatcher.matchBiopsies(patient.patientId(),
                sequencedBiopsies, clinicalBiopsies);
        ThreadContext.put("cpctHospitalCode", "default");
        return new Patient(patientInfo, sequencedBiopsies, matchedTreatments, matchedResponses, matchedBiopsies);
    }
}
