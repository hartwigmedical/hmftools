package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientData;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.data.TumorMarkerData;
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
    private final LimsSampleReader limsSampleReader;
    @NotNull
    private final PreTreatmentReader preTreatmentReader;
    @NotNull
    private final BiopsyTreatmentReader biopsyTreatmentReader;

    public PatientReader(@NotNull final CpctEcrfModel model, @NotNull final TreatmentCurator treatmentCurator,
            @NotNull final TumorLocationCurator tumorLocationCurator, @NotNull final Lims lims) {
        cpctPatientReader = new CpctPatientReader(model, tumorLocationCurator);
        limsSampleReader = new LimsSampleReader(lims);
        preTreatmentReader = new PreTreatmentReader(treatmentCurator);
        biopsyTreatmentReader = new BiopsyTreatmentReader(treatmentCurator);
    }

    @NotNull
    public Patient read(@NotNull final EcrfPatient ecrfPatient, @NotNull final List<String> tumorSamplesForPatient) throws IOException {
        LOGGER.info("Analyzing patient " + ecrfPatient.patientId() + " with samples: " + tumorSamplesForPatient);

        final List<SampleData> sequencedBiopsies = limsSampleReader.read(tumorSamplesForPatient);
        final PatientData patientData = cpctPatientReader.read(ecrfPatient);
        final PreTreatmentData preTreatmentData = preTreatmentReader.read(ecrfPatient);
        final List<BiopsyData> clinicalBiopsies = BiopsyReader.read(ecrfPatient);
        final List<BiopsyTreatmentData> treatments = biopsyTreatmentReader.read(ecrfPatient);
        final List<BiopsyTreatmentResponseData> treatmentResponses = BiopsyTreatmentResponseReader.read(ecrfPatient);
        final List<TumorMarkerData> tumorMarkers = TumorMarkerReader.read(ecrfPatient);

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

        return new Patient(patientData,
                preTreatmentData,
                sequencedBiopsies,
                matchedBiopsies.values(),
                matchedTreatments.values(), matchedResponses.values(), tumorMarkers,
                findings);
    }
}
