package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyTreatmentBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyTreatmentResponseBuilder;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentReader.FORM_TREATMENT;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentResponseReader.FIELD_ASSESSMENT_DATE;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentResponseReader.FIELD_MEASUREMENT_DONE;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentResponseReader.FIELD_RESPONSE;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentResponseReader.FIELD_RESPONSE_DATE;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentResponseReader.FORM_TUMOR_MEASUREMENT;
import static com.hartwig.hmftools.patientdb.validators.PatientValidator.fields;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.ImmutableDrugData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class TreatmentResponseValidationTest {
    private final String CPCT_ID = "CPCT01020000";
    private final static LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private final static LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private final static LocalDate MAR2015 = LocalDate.parse("2015-03-01");

    private final static BiopsyTreatmentData TREATMENT_JAN_MAR =
            biopsyTreatmentBuilder().addDrugs(drugWithStartAndEndDate(JAN2015, MAR2015)).build();
    private final static BiopsyTreatmentData TREATMENT_JAN_ONGOING =
            biopsyTreatmentBuilder().addDrugs(drugWithStartAndEndDate(JAN2015, null)).build();
    private final static BiopsyTreatmentData TREATMENT_JAN_JAN =
            biopsyTreatmentBuilder().addDrugs(drugWithStartAndEndDate(JAN2015, JAN2015)).build();

    private final static BiopsyTreatmentResponseData RESPONSE_JAN2015 =
            biopsyTreatmentResponseBuilder().measurementDone("Yes").response("PR").responseDate(JAN2015).build();
    private final static BiopsyTreatmentResponseData RESPONSE_FEB2015 =
            biopsyTreatmentResponseBuilder().measurementDone("Yes").response("PR").responseDate(FEB2015).build();
    private final static BiopsyTreatmentResponseData RESPONSE_NULL = biopsyTreatmentResponseBuilder().build();
    private final static BiopsyTreatmentResponseData RESPONSE_ONLY = biopsyTreatmentResponseBuilder().response("PR").build();
    private final static BiopsyTreatmentResponseData RESPONSE_MISSING_DATA =
            biopsyTreatmentResponseBuilder().measurementDone("Yes").build();
    private final static BiopsyTreatmentResponseData RESPONSE_MEASUREMENT_NO_WITH_DATA =
            biopsyTreatmentResponseBuilder().measurementDone("No").response("PR").responseDate(JAN2015).build();
    private final static BiopsyTreatmentResponseData RESPONSE_MEASUREMENT_NO_WITH_VALID_DATA =
            biopsyTreatmentResponseBuilder().measurementDone("No").response("ND").responseDate(JAN2015).build();

    @Test
    public void reportsMeasurementDoneNull() {
        final List<ValidationFinding> findings = validateNonFirstResponse(CPCT_ID, RESPONSE_NULL);
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_MEASUREMENT_DONE));
    }

    @Test
    public void reportsOnlyResponseFilledIn() {
        final List<ValidationFinding> findings = validateNonFirstResponse(CPCT_ID, RESPONSE_ONLY);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_MEASUREMENT_DONE));
        assertTrue(findingsFields.contains(fields(FIELD_ASSESSMENT_DATE, FIELD_RESPONSE_DATE)));
    }

    @Test
    public void reportsMeasurementDoneMissingData() {
        final List<ValidationFinding> findings = validateNonFirstResponse(CPCT_ID, RESPONSE_MISSING_DATA);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_RESPONSE));
        assertTrue(findingsFields.contains(fields(FIELD_ASSESSMENT_DATE, FIELD_RESPONSE_DATE)));
    }

    @Test
    public void acceptMissingDataForFirstResponse() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponse(CPCT_ID, RESPONSE_MISSING_DATA, true);
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertFalse(findingsFields.contains(FIELD_RESPONSE));
        assertTrue(findingsFields.contains(fields(FIELD_ASSESSMENT_DATE, FIELD_RESPONSE_DATE)));
    }

    @Test
    public void reportsMeasurementDoneNoWithData() {
        final List<ValidationFinding> findings = validateNonFirstResponse(CPCT_ID, RESPONSE_MEASUREMENT_NO_WITH_DATA);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_MEASUREMENT_DONE));
    }

    @Test
    public void ignoresMeasurementDoneNoWithValidData() {
        final List<ValidationFinding> findings = validateNonFirstResponse(CPCT_ID, RESPONSE_MEASUREMENT_NO_WITH_VALID_DATA);
        assertTrue(findings.isEmpty());
    }

    @Test
    public void reportsFirstMeasurementAfterTreatmentStart() {
        BiopsyTreatmentResponseData matchedResponseFeb2015 =
                ImmutableBiopsyTreatmentResponseData.builder().from(RESPONSE_FEB2015).treatmentId(TREATMENT_JAN_MAR.id()).build();
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponses(CPCT_ID,
                Lists.newArrayList(TREATMENT_JAN_MAR),
                Lists.newArrayList(matchedResponseFeb2015));
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(fields(FORM_TREATMENT, FORM_TUMOR_MEASUREMENT)));
    }

    @Test
    public void reportsMissingTreatmentData() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatmentResponses(CPCT_ID, Lists.newArrayList(), Lists.newArrayList(RESPONSE_JAN2015));
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FORM_TREATMENT));
    }

    @Test
    public void reportsMissingResponseForTreatmentMoreThan16Weeks() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatmentResponses(CPCT_ID, Lists.newArrayList(TREATMENT_JAN_ONGOING), Lists.newArrayList());
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FORM_TUMOR_MEASUREMENT));
    }

    @Test
    public void doesNotReportMissingResponseForTreatmentEndedImmediately() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatmentResponses(CPCT_ID, Lists.newArrayList(TREATMENT_JAN_JAN), Lists.newArrayList());
        assertEquals(0, findings.size());
    }

    @NotNull
    private static DrugData drugWithStartAndEndDate(@Nullable LocalDate startDate, @Nullable LocalDate endDate) {
        return ImmutableDrugData.of("anything", startDate, endDate, null, Lists.newArrayList());
    }

    @NotNull
    private static List<ValidationFinding> validateNonFirstResponse(@NotNull String patientId,
            @NotNull BiopsyTreatmentResponseData response) {
        return PatientValidator.validateTreatmentResponse(patientId, response, false);
    }
}
