package com.hartwig.hmftools.patientdb.clinical.readers.cpct;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.consents.ConsentConfig;
import com.hartwig.hmftools.patientdb.clinical.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BaselineData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.PreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.RanoMeasurementData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.TumorMarkerData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.clinical.matchers.BiopsyMatcher;
import com.hartwig.hmftools.patientdb.clinical.matchers.MatchResult;
import com.hartwig.hmftools.patientdb.clinical.matchers.TreatmentMatcher;
import com.hartwig.hmftools.patientdb.clinical.matchers.TreatmentResponseMatcher;
import com.hartwig.hmftools.patientdb.clinical.readers.EcrfPatientReader;

import org.jetbrains.annotations.NotNull;

public class CpctPatientReader implements EcrfPatientReader {

    @NotNull
    private final BaselineReader baselineReader;
    @NotNull
    private final PreTreatmentReader preTreatmentReader;
    @NotNull
    private final BiopsyReader biopsyReader;
    @NotNull
    private final BiopsyTreatmentReader biopsyTreatmentReader;

    public CpctPatientReader(@NotNull PrimaryTumorCurator primaryTumorCurator, @NotNull Map<Integer, String> hospitals,
            @NotNull BiopsySiteCurator biopsySiteCurator, @NotNull TreatmentCurator treatmentCurator) {
        this.baselineReader = new BaselineReader(primaryTumorCurator, hospitals);
        this.preTreatmentReader = new PreTreatmentReader(treatmentCurator);
        this.biopsyReader = new BiopsyReader(biopsySiteCurator);
        this.biopsyTreatmentReader = new BiopsyTreatmentReader(treatmentCurator);
    }

    @NotNull
    @Override
    public Patient read(@NotNull EcrfPatient ecrfPatient, @NotNull List<SampleData> sequencedSamples,
            @NotNull Map<String, ConsentConfig> consentConfigMap) throws IOException {
        BaselineData baselineData = baselineReader.read(ecrfPatient, consentConfigMap);
        PreTreatmentData preTreatmentData = preTreatmentReader.read(ecrfPatient);
        List<BiopsyData> clinicalBiopsies = biopsyReader.read(ecrfPatient, baselineData.curatedPrimaryTumor());
        List<BiopsyTreatmentData> treatments = biopsyTreatmentReader.read(ecrfPatient);
        List<BiopsyTreatmentResponseData> treatmentResponses = BiopsyTreatmentResponseReader.read(ecrfPatient);
        List<TumorMarkerData> tumorMarkers = TumorMarkerReader.read(ecrfPatient);
        List<RanoMeasurementData> ranoMeasurements = RanoMeasurementReader.read(ecrfPatient);

        MatchResult<BiopsyData> matchedBiopsies =
                BiopsyMatcher.matchBiopsiesToTumorSamples(ecrfPatient.patientId(), sequencedSamples, clinicalBiopsies);
        MatchResult<BiopsyTreatmentData> matchedTreatments =
                TreatmentMatcher.matchTreatmentsToBiopsies(ecrfPatient.patientId(), withSampleMatchOnly(matchedBiopsies), treatments);

        // We also match responses to unmatched treatments. Not sure that is optimal. See also DEV-477.
        MatchResult<BiopsyTreatmentResponseData> matchedResponses =
                TreatmentResponseMatcher.matchTreatmentResponsesToTreatments(ecrfPatient.patientId(),
                        matchedTreatments.values(),
                        treatmentResponses);

        List<ValidationFinding> findings = Lists.newArrayList();
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
                ranoMeasurements,
                findings);
    }

    @NotNull
    private static List<BiopsyData> withSampleMatchOnly(@NotNull MatchResult<BiopsyData> biopsies) {
        List<BiopsyData> biopsiesWithMatchedSample = Lists.newArrayList();
        for (BiopsyData biopsy : biopsies.values()) {
            if (biopsy.sampleId() != null) {
                biopsiesWithMatchedSample.add(biopsy);
            }
        }
        return biopsiesWithMatchedSample;
    }
}
