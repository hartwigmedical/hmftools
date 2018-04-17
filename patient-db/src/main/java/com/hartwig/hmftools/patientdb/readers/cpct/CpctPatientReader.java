package com.hartwig.hmftools.patientdb.readers.cpct;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.curators.BiopsySiteCurator;
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
import com.hartwig.hmftools.patientdb.readers.PatientReader;

import org.jetbrains.annotations.NotNull;

public class CpctPatientReader implements PatientReader {
    @NotNull
    private final BaselineReader baselineReader;
    @NotNull
    private final PreTreatmentReader preTreatmentReader;
    @NotNull
    private final BiopsyReader biopsyReader;
    @NotNull
    private final BiopsyTreatmentReader biopsyTreatmentReader;

    public CpctPatientReader(@NotNull TumorLocationCurator tumorLocationCurator, @NotNull Map<Integer, String> hospitals,
            @NotNull BiopsySiteCurator biopsySiteCurator, @NotNull TreatmentCurator treatmentCurator) {
        this.baselineReader = new BaselineReader(tumorLocationCurator, hospitals);
        this.preTreatmentReader = new PreTreatmentReader(treatmentCurator);
        this.biopsyReader = new BiopsyReader(biopsySiteCurator);
        this.biopsyTreatmentReader = new BiopsyTreatmentReader(treatmentCurator);
    }

    @NotNull
    @Override
    public Patient read(@NotNull final EcrfPatient ecrfPatient, @NotNull final List<SampleData> sequencedBiopsies) {
        final BaselineData baselineData = baselineReader.read(ecrfPatient);
        final PreTreatmentData preTreatmentData = preTreatmentReader.read(ecrfPatient);
        final List<BiopsyData> clinicalBiopsies = biopsyReader.read(ecrfPatient, baselineData.curatedTumorLocation());
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

        return new Patient(ecrfPatient.patientId(),
                baselineData,
                preTreatmentData, sequencedBiopsies,
                matchedBiopsies.values(),
                matchedTreatments.values(),
                matchedResponses.values(),
                tumorMarkers,
                findings);
    }
}
