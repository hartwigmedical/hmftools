package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentReader.FORM_TREATMENT;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentResponseReader.FIELD_ASSESSMENT_DATE;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentResponseReader.FIELD_MEASUREMENT_YN;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentResponseReader.FIELD_RESPONSE;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentResponseReader.FORM_TUMOR_MEASUREMENT;
import static com.hartwig.hmftools.patientdb.validators.PatientValidator.fields;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;

import org.junit.Test;

public class TreatmentResponseValidationTest {
    private final String CPCT_ID = "CPCT01020000";
    private final static LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private final static LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private final static LocalDate MAR2015 = LocalDate.parse("2015-03-01");

    private final static BiopsyTreatmentData TREATMENT_JAN_MAR = new BiopsyTreatmentData("Yes", JAN2015, MAR2015, Lists.newArrayList());
    private final static BiopsyTreatmentData TREATMENT_JAN_ONGOING = new BiopsyTreatmentData("Yes", JAN2015, null, Lists.newArrayList());
    private final static BiopsyTreatmentData TREATMENT_JAN_JAN = new BiopsyTreatmentData("Yes", JAN2015, JAN2015, Lists.newArrayList());
    private final static BiopsyTreatmentResponseData RESPONSE_JAN2015 = new BiopsyTreatmentResponseData(JAN2015, JAN2015, "NE", "Yes");
    private final static BiopsyTreatmentResponseData RESPONSE_FEB2015 = new BiopsyTreatmentResponseData(FEB2015, FEB2015, "NE", "Yes");
    private final static BiopsyTreatmentResponseData RESPONSE_NULL = new BiopsyTreatmentResponseData(null, null, null, null);
    private final static BiopsyTreatmentResponseData RESPONSE_ONLY = new BiopsyTreatmentResponseData(null, null, "NE", null);
    private final static BiopsyTreatmentResponseData RESPONSE_MISSING_DATA = new BiopsyTreatmentResponseData(null, null, null, "yes");
    private final static BiopsyTreatmentResponseData RESPONSE_MEASUREMENT_NO_WITH_DATA =
            new BiopsyTreatmentResponseData(JAN2015, JAN2015, "NE", "no");

    @Test
    public void reportsMeasurementDoneNull() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponse(CPCT_ID, RESPONSE_NULL);
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_MEASUREMENT_YN));
    }

    @Test
    public void reportsOnlyResponseFilledIn() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponse(CPCT_ID, RESPONSE_ONLY);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_MEASUREMENT_YN));
        assertTrue(findingsFields.contains(FIELD_ASSESSMENT_DATE));
    }

    @Test
    public void reportsMeasurementDoneMissingData() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponse(CPCT_ID, RESPONSE_MISSING_DATA);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_RESPONSE));
        assertTrue(findingsFields.contains(FIELD_ASSESSMENT_DATE));
    }

    @Test
    public void reportsMeasurementDoneNoWithData() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponse(CPCT_ID, RESPONSE_MEASUREMENT_NO_WITH_DATA);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_MEASUREMENT_YN));
    }

    @Test
    public void reportsFirstMeasurementAfterTreatmentStart() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponses(CPCT_ID, Lists.newArrayList(TREATMENT_JAN_MAR),
                Lists.newArrayList(RESPONSE_FEB2015));
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
}
