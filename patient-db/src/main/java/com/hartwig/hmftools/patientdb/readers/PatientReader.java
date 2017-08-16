package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientData;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.matchers.BiopsyMatcher;
import com.hartwig.hmftools.patientdb.matchers.MatchResult;
import com.hartwig.hmftools.patientdb.matchers.TreatmentMatcher;
import com.hartwig.hmftools.patientdb.matchers.TreatmentResponseMatcher;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PatientReader {
    private static final Logger LOGGER = LogManager.getLogger(PatientReader.class);

    @NotNull
    private final CpctPatientReader cpctPatientReader;
    @NotNull
    private final SampleReader sampleReader;
    @NotNull
    private final BiopsyTreatmentReader biopsyTreatmentReader;

    public PatientReader(@NotNull final CpctEcrfModel model, @NotNull final Map<String, String> treatmentToTypeMappings,
            @NotNull final String limsCsv, @NotNull final String limsOldCsv, @NotNull final String umcuCsv)
            throws IOException, HartwigException {
        cpctPatientReader = new CpctPatientReader(model);
        sampleReader = new SampleReader(limsCsv, limsOldCsv, umcuCsv);
        biopsyTreatmentReader = new BiopsyTreatmentReader(treatmentToTypeMappings);
    }

    @NotNull
    public Patient read(@NotNull final EcrfPatient ecrfPatient, @NotNull final List<String> tumorSamplesForPatient)
            throws IOException, HartwigException {
        LOGGER.info("Reading patient " + ecrfPatient.patientId() + ": Samples: " + tumorSamplesForPatient);
        final List<SampleData> sequencedBiopsies = sampleReader.read(tumorSamplesForPatient);
        final PatientData patientData = cpctPatientReader.read(ecrfPatient);
        final List<BiopsyData> clinicalBiopsies = BiopsyReader.read(ecrfPatient);
        final List<BiopsyTreatmentData> treatments = biopsyTreatmentReader.read(ecrfPatient);
        final List<BiopsyTreatmentResponseData> treatmentResponses = BiopsyTreatmentResponseReader.read(ecrfPatient);
        final MatchResult<BiopsyData> matchedBiopsies =
                BiopsyMatcher.matchBiopsiesToTumorSamples(ecrfPatient.patientId(), sequencedBiopsies, clinicalBiopsies);
        final MatchResult<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies(ecrfPatient.patientId(), clinicalBiopsies, treatments);
        final MatchResult<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments(ecrfPatient.patientId(), treatments, treatmentResponses);
        final List<ValidationFinding> findings = Lists.newArrayList();
        findings.addAll(matchedBiopsies.findings());
        findings.addAll(matchedTreatments.findings());
        findings.addAll(matchedResponses.findings());
        return new Patient(patientData, sequencedBiopsies, matchedBiopsies.values(), matchedTreatments.values(), matchedResponses.values(),
                findings);
    }
}
