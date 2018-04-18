package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.baselineBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyTreatmentBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyTreatmentResponseBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.drugBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.CuratedTreatment;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTreatment;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocation;
import com.hartwig.hmftools.patientdb.data.ImmutableDrugData;
import com.hartwig.hmftools.patientdb.data.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class PatientValidatorTest {

    private static final String PATIENT_IDENTIFIER = "CPCT01020000";

    private static final LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private static final LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private static final LocalDate MAR2015 = LocalDate.parse("2015-03-01");

    private static final String HOSPITAL = "Test Hospital";
    private static final BaselineData EMPTY_BASELINE = baselineBuilder().build();
    private static final BaselineData BASELINE_DATA_MISSING_LOCATION_MAPPING = baselineBuilder().hospital(HOSPITAL)
            .curatedTumorLocation(ImmutableCuratedTumorLocation.of(null, null, "some_location"))
            .build();

    private static final BiopsyData BIOPSY_NULL = biopsyBuilder().date(null).build();
    private static final BiopsyData BIOPSY_FEB1 = biopsyBuilder().date(FEB2015).site("1").location("").build();
    private static final BiopsyData BIOPSY_FEB2 = biopsyBuilder().date(FEB2015).site("2").location("").build();

    private static final DrugData DRUG_NULL = create(null, null, null);
    private static final DrugData DRUG_WRONG = create(null, FEB2015, JAN2015);
    private static final DrugData DRUG_JAN_JAN = create("Drug1", JAN2015, JAN2015);
    private static final DrugData DRUG_JAN_ONGOING = create("Drug1", JAN2015, null);
    private static final DrugData DRUG_JAN_FEB = create("Drug1", JAN2015, FEB2015);
    private static final DrugData DRUG_FEB_ONGOING = create("Drug1", FEB2015, null);
    private static final DrugData DRUG_JAN_MAR = create("Drug1", JAN2015, MAR2015);

    private static final DrugData DRUG_WITH_PARTIAL_CURATED_ENTRY = ImmutableDrugData.of("Drug1 Drug2 Drug3",
            JAN2015,
            JAN2015,
            null,
            Lists.newArrayList(ImmutableCuratedTreatment.of("Drug1", "Type1", "Drug1")));
    private static final DrugData DRUG_MISSING_CURATED_ENTRY = ImmutableDrugData.of("Drug1", JAN2015, JAN2015, null, Lists.newArrayList());

    private static final BiopsyTreatmentData TREATMENT_RADIO_THERAPY_NULL =
            biopsyTreatmentBuilder().treatmentGiven("No").radiotherapyGiven(null).build();
    private static final BiopsyTreatmentData TREATMENT_GIVEN_NULL =
            biopsyTreatmentBuilder().treatmentGiven(null).radiotherapyGiven("No").build();
    private static final BiopsyTreatmentData TREATMENT_GIVEN_EMPTY = biopsyTreatmentBuilder().treatmentGiven("Yes").build();
    private static final BiopsyTreatmentData TREATMENT_NOT_GIVEN_DATA =
            biopsyTreatmentBuilder().treatmentGiven("No").addDrugs(DRUG_JAN_FEB).radiotherapyGiven("Yes").build();
    private static final BiopsyTreatmentData TREATMENT_GIVEN_GIBBERISH =
            biopsyTreatmentBuilder().treatmentGiven("mmm").radiotherapyGiven("Yes").build();
    private static final BiopsyTreatmentData TREATMENT_WRONG_DRUG_DATA =
            biopsyTreatmentBuilder().treatmentGiven("Yes").radiotherapyGiven("Yes").addDrugs(DRUG_NULL, DRUG_WRONG).build();

    private static final BiopsyTreatmentData TREATMENT_JAN_JAN =
            biopsyTreatmentBuilder().treatmentGiven("Yes").radiotherapyGiven("Yes").addDrugs(DRUG_JAN_JAN).build();
    private static final BiopsyTreatmentData TREATMENT_JAN_FEB =
            biopsyTreatmentBuilder().treatmentGiven("Yes").radiotherapyGiven("Yes").addDrugs(DRUG_JAN_FEB).build();
    private static final BiopsyTreatmentData TREATMENT_JAN_MAR =
            biopsyTreatmentBuilder().treatmentGiven("Yes").radiotherapyGiven("Yes").addDrugs(DRUG_JAN_MAR).build();
    private static final BiopsyTreatmentData TREATMENT_JAN_ONGOING =
            biopsyTreatmentBuilder().treatmentGiven("Yes").radiotherapyGiven("Yes").addDrugs(DRUG_JAN_ONGOING).build();
    private static final BiopsyTreatmentData TREATMENT_FEB_ONGOING =
            biopsyTreatmentBuilder().treatmentGiven("Yes").radiotherapyGiven("Yes").addDrugs(DRUG_FEB_ONGOING).build();

    private static final BiopsyTreatmentResponseData RESPONSE_JAN2015 =
            biopsyTreatmentResponseBuilder().measurementDone("Yes").response("PR").responseDate(JAN2015).build();
    private static final BiopsyTreatmentResponseData RESPONSE_FEB2015 =
            biopsyTreatmentResponseBuilder().measurementDone("Yes").response("PR").responseDate(FEB2015).build();
    private static final BiopsyTreatmentResponseData RESPONSE_NULL = biopsyTreatmentResponseBuilder().build();
    private static final BiopsyTreatmentResponseData RESPONSE_ONLY = biopsyTreatmentResponseBuilder().response("PR").build();
    private static final BiopsyTreatmentResponseData RESPONSE_MISSING_DATA =
            biopsyTreatmentResponseBuilder().measurementDone("Yes").build();
    private static final BiopsyTreatmentResponseData RESPONSE_MEASUREMENT_NO_WITH_DATA =
            biopsyTreatmentResponseBuilder().measurementDone("No").response("PR").responseDate(JAN2015).build();
    private static final BiopsyTreatmentResponseData RESPONSE_MEASUREMENT_NO_WITH_VALID_DATA =
            biopsyTreatmentResponseBuilder().measurementDone("No").response("ND").responseDate(JAN2015).build();

    @Test
    public void reportsMissingBaselineFields() {
        final List<ValidationFinding> findings = PatientValidator.validateBaselineData(PATIENT_IDENTIFIER, EMPTY_BASELINE);
        assertEquals(6, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void reportMissingPreTreatment() {
        PreTreatmentData emptyData = ImmutablePreTreatmentData.builder().formStatus(FormStatus.undefined()).build();

        assertEquals(2, PatientValidator.validatePreTreatmentData(PATIENT_IDENTIFIER, emptyData).size());

        PreTreatmentData actualData = ImmutablePreTreatmentData.builder()
                .radiotherapyGiven("Yes")
                .treatmentGiven("No")
                .formStatus(FormStatus.undefined())
                .build();
        assertEquals(0, PatientValidator.validatePreTreatmentData(PATIENT_IDENTIFIER, actualData).size());

        PreTreatmentData onlyRadioTherapyPresent =
                ImmutablePreTreatmentData.builder().radiotherapyGiven("Yes").formStatus(FormStatus.undefined()).build();
        assertEquals(1, PatientValidator.validatePreTreatmentData(PATIENT_IDENTIFIER, onlyRadioTherapyPresent).size());

        PreTreatmentData onlyTreatmentPresent =
                ImmutablePreTreatmentData.builder().treatmentGiven("Yes").formStatus(FormStatus.undefined()).build();
        assertEquals(1, PatientValidator.validatePreTreatmentData(PATIENT_IDENTIFIER, onlyTreatmentPresent).size());
    }

    @Test
    public void reportsFailureToCuratePrimaryTumorLocation() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTumorLocationCuration(PATIENT_IDENTIFIER, BASELINE_DATA_MISSING_LOCATION_MAPPING);
        assertEquals(1, findings.size());
        assertEquals("tumorLocationCuration", findings.get(0).level());
    }

    @Test
    public void reportsMissingBiopsyFields() {
        final List<ValidationFinding> findings = PatientValidator.validateBiopsyData(PATIENT_IDENTIFIER, BIOPSY_NULL);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void reportsAllBiopsyFieldsEmpty() {
        final List<ValidationFinding> findings =
                PatientValidator.validateBiopsies(PATIENT_IDENTIFIER, Lists.newArrayList(BIOPSY_NULL, BIOPSY_FEB1, BIOPSY_FEB2));
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void reportsBiopsyBeforeInformedConsent() {
        final List<ValidationFinding> findings = PatientValidator.validateInformedConsentDate(PATIENT_IDENTIFIER,
                baselineBuilder().informedConsentDate(MAR2015).build(),
                Lists.newArrayList(BIOPSY_FEB1));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsMissingDrugData() {
        final List<ValidationFinding> findings = PatientValidator.validateDrugData(PATIENT_IDENTIFIER, DRUG_NULL, FormStatus.undefined());
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void reportsIncorrectDrugData() {
        final List<ValidationFinding> findings = PatientValidator.validateDrugData(PATIENT_IDENTIFIER, DRUG_WRONG, FormStatus.undefined());
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void reportsMissingRadioTherapy() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(PATIENT_IDENTIFIER, Lists.newArrayList(TREATMENT_RADIO_THERAPY_NULL));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsMissingTreatmentGiven() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(PATIENT_IDENTIFIER, Lists.newArrayList(TREATMENT_GIVEN_NULL));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsMissingTreatmentData() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(PATIENT_IDENTIFIER, Lists.newArrayList(TREATMENT_GIVEN_EMPTY));
        assertEquals(2, findings.size());
    }

    @Test
    public void reportsTreatmentGivenNoWithData() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(PATIENT_IDENTIFIER, Lists.newArrayList(TREATMENT_NOT_GIVEN_DATA));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsTreatmentGivenGibberish() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(PATIENT_IDENTIFIER, Lists.newArrayList(TREATMENT_GIVEN_GIBBERISH));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsDrugFindingsForTreatment() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(PATIENT_IDENTIFIER, Lists.newArrayList(TREATMENT_WRONG_DRUG_DATA));
        assertEquals(4, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void reportsWrongTreatmentTimeline() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatments(PATIENT_IDENTIFIER,
                Lists.newArrayList(TREATMENT_JAN_MAR, TREATMENT_JAN_JAN, TREATMENT_JAN_FEB));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsWrongTreatmentTimelineOngoing() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(PATIENT_IDENTIFIER, Lists.newArrayList(TREATMENT_JAN_FEB, TREATMENT_JAN_ONGOING));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsTwoOngoingTreatments() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(PATIENT_IDENTIFIER, Lists.newArrayList(TREATMENT_JAN_ONGOING, TREATMENT_FEB_ONGOING));
        assertEquals(3, findings.size());
    }

    @Test
    public void reportsMissingCuratedTreatment() {
        String curationName = "testTreatmentCuration";
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentCuration(PATIENT_IDENTIFIER,
                curationName,
                Lists.newArrayList(ImmutableBiopsyTreatmentData.of("Yes",
                        "Yes",
                        Lists.newArrayList(DRUG_MISSING_CURATED_ENTRY),
                        FormStatus.undefined())));
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::level).collect(Collectors.toList());
        assertTrue(findingsFields.get(0).equals(curationName));
    }

    @Test
    public void reportsPartiallyCuratedTreatment() {
        String curationName = "testTreatmentCuration";
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentCuration(PATIENT_IDENTIFIER,
                curationName,
                Lists.newArrayList(ImmutableBiopsyTreatmentData.of("Yes",
                        "Yes",
                        Lists.newArrayList(DRUG_WITH_PARTIAL_CURATED_ENTRY),
                        FormStatus.undefined())));
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::level).collect(Collectors.toList());
        assertTrue(findingsFields.get(0).equals(curationName));
    }

    @Test
    public void doesNotReportCorrectDeathTimeline() {
        final List<ValidationFinding> findings = PatientValidator.validateDeathDate(PATIENT_IDENTIFIER,
                baselineBuilder().deathDate(MAR2015).build(),
                Lists.newArrayList(TREATMENT_JAN_JAN, TREATMENT_JAN_FEB));
        assertEquals(0, findings.size());
    }

    @Test
    public void reportsDeathDateBeforeEndOfTreatment() {
        final List<ValidationFinding> findings = PatientValidator.validateDeathDate(PATIENT_IDENTIFIER,
                baselineBuilder().deathDate(MAR2015).build(),
                Lists.newArrayList(TREATMENT_JAN_ONGOING, TREATMENT_JAN_FEB));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsMeasurementDoneNull() {
        final List<ValidationFinding> findings = validateNonFirstResponse(PATIENT_IDENTIFIER, RESPONSE_NULL);
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsOnlyResponseFilledIn() {
        final List<ValidationFinding> findings = validateNonFirstResponse(PATIENT_IDENTIFIER, RESPONSE_ONLY);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void reportsMeasurementDoneMissingData() {
        final List<ValidationFinding> findings = validateNonFirstResponse(PATIENT_IDENTIFIER, RESPONSE_MISSING_DATA);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void acceptMissingDataForFirstResponse() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatmentResponse(PATIENT_IDENTIFIER, RESPONSE_MISSING_DATA, true);
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void reportsMeasurementDoneNoWithData() {
        final List<ValidationFinding> findings = validateNonFirstResponse(PATIENT_IDENTIFIER, RESPONSE_MEASUREMENT_NO_WITH_DATA);
        assertEquals(2, findings.size());
    }

    @Test
    public void ignoresMeasurementDoneNoWithValidData() {
        final List<ValidationFinding> findings = validateNonFirstResponse(PATIENT_IDENTIFIER, RESPONSE_MEASUREMENT_NO_WITH_VALID_DATA);
        assertTrue(findings.isEmpty());
    }

    @Test
    public void reportsFirstMeasurementAfterTreatmentStart() {
        BiopsyTreatmentResponseData matchedResponseFeb2015 =
                ImmutableBiopsyTreatmentResponseData.builder().from(RESPONSE_FEB2015).treatmentId(TREATMENT_JAN_MAR.id()).build();
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponses(PATIENT_IDENTIFIER,
                Lists.newArrayList(TREATMENT_JAN_MAR),
                Lists.newArrayList(matchedResponseFeb2015));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsMissingTreatmentForResponseData() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatmentResponses(PATIENT_IDENTIFIER, Lists.newArrayList(), Lists.newArrayList(RESPONSE_JAN2015));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsMissingResponseForTreatmentMoreThan16Weeks() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponses(PATIENT_IDENTIFIER,
                Lists.newArrayList(TREATMENT_JAN_ONGOING),
                Lists.newArrayList());
        assertEquals(1, findings.size());
    }

    @Test
    public void doesNotReportMissingResponseForTreatmentEndedImmediately() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponses(PATIENT_IDENTIFIER,
                Lists.newArrayList(TREATMENT_JAN_JAN),
                Lists.newArrayList());
        assertEquals(0, findings.size());
    }

    @NotNull
    private static DrugData create(@Nullable String name, @Nullable LocalDate startDate, @Nullable LocalDate endDate) {
        List<CuratedTreatment> curation =
                name != null ? Lists.newArrayList(ImmutableCuratedTreatment.of(name, "Type1", name)) : Lists.newArrayList();
        return drugBuilder().name(name).startDate(startDate).endDate(endDate).addAllCuratedTreatments(curation).build();
    }

    @NotNull
    private static List<ValidationFinding> validateNonFirstResponse(@NotNull String patientId,
            @NotNull BiopsyTreatmentResponseData response) {
        return PatientValidator.validateTreatmentResponse(patientId, response, false);
    }
}