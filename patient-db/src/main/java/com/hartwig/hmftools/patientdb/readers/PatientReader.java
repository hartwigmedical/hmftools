package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.Patient;
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
    private final PreTreatmentReader preTreatmentReader;
    @NotNull
    private final BiopsyTreatmentReader biopsyTreatmentReader;

    public PatientReader(@NotNull final CpctEcrfModel model, @NotNull final TreatmentCurator treatmentCurator,
            @NotNull final TumorLocationCurator tumorLocationCurator) {
        this.cpctPatientReader = new CpctPatientReader(model, tumorLocationCurator);
        this.preTreatmentReader = new PreTreatmentReader(treatmentCurator);
        this.biopsyTreatmentReader = new BiopsyTreatmentReader(treatmentCurator);
    }

    @NotNull
    public Patient read(@NotNull final EcrfPatient ecrfPatient, @NotNull final List<SampleData> sequencedSamples) throws IOException {
        LOGGER.info("Reading patient " + ecrfPatient.patientId() + " with samples: " + sequencedSamples);

        final BaselineData baselineData = cpctPatientReader.read(ecrfPatient);
        final PreTreatmentData preTreatmentData = preTreatmentReader.read(ecrfPatient);
        final List<BiopsyData> clinicalBiopsies = BiopsyReader.read(ecrfPatient);
        final List<BiopsyTreatmentData> treatments = biopsyTreatmentReader.read(ecrfPatient);
        final List<BiopsyTreatmentResponseData> treatmentResponses = BiopsyTreatmentResponseReader.read(ecrfPatient);
        final List<TumorMarkerData> tumorMarkers = TumorMarkerReader.read(ecrfPatient);

        final MatchResult<BiopsyData> matchedBiopsies =
                BiopsyMatcher.matchBiopsiesToTumorSamples(ecrfPatient.patientId(), sequencedSamples, clinicalBiopsies);
        final MatchResult<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies(ecrfPatient.patientId(), clinicalBiopsies, treatments);
        final MatchResult<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments(ecrfPatient.patientId(), treatments, treatmentResponses);

        final List<ValidationFinding> findings = Lists.newArrayList();
        findings.addAll(matchedBiopsies.findings());
        findings.addAll(matchedTreatments.findings());
        findings.addAll(matchedResponses.findings());

        return new Patient(ecrfPatient.patientId(),
                baselineData,
                preTreatmentData,
                sequencedSamples,
                matchedBiopsies.values(),
                matchedTreatments.values(),
                matchedResponses.values(),
                tumorMarkers,
                findings);
    }
}
