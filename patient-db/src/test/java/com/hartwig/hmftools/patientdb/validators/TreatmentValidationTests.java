package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentReader.FIELD_DRUG;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentReader.FIELD_DRUG_END;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentReader.FIELD_DRUG_START;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentReader.FIELD_TREATMENT_GIVEN;
import static com.hartwig.hmftools.patientdb.readers.BiopsyTreatmentReader.FORM_TREATMENT;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_DEATH_DATE;
import static com.hartwig.hmftools.patientdb.validators.PatientValidator.fields;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentDrugData;

import org.junit.Test;

public class TreatmentValidationTests {
    private final String CPCT_ID = "CPCT01020000";
    private final static LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private final static LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private final static LocalDate MAR2015 = LocalDate.parse("2015-03-01");

    private final static BiopsyTreatmentDrugData DRUG_NULL = new BiopsyTreatmentDrugData(null, null, null, null);
    private final static BiopsyTreatmentDrugData DRUG_WRONG = new BiopsyTreatmentDrugData(null, null, FEB2015, JAN2015);
    private final static BiopsyTreatmentDrugData DRUG_JAN_FEB = new BiopsyTreatmentDrugData("Drug1", "Type1", JAN2015, FEB2015);
    private final static BiopsyTreatmentDrugData DRUG_JAN_MAR = new BiopsyTreatmentDrugData("Drug1", "Type1", JAN2015, MAR2015);
    private final static BiopsyTreatmentDrugData DRUG_JAN_JAN = new BiopsyTreatmentDrugData("Drug1", "Type1", JAN2015, JAN2015);

    private final static BiopsyTreatmentData TREATMENT_GIVEN_NULL = new BiopsyTreatmentData(null, null, null, Lists.newArrayList());
    private final static BiopsyTreatmentData TREATMENT_GIVEN_EMPTY = new BiopsyTreatmentData("Yes", null, null, Lists.newArrayList());
    private final static BiopsyTreatmentData TREATMENT_NOTGIVEN_DATA =
            new BiopsyTreatmentData("No", JAN2015, FEB2015, Lists.newArrayList(DRUG_JAN_FEB));
    private final static BiopsyTreatmentData TREATMENT_GIVEN_GIBBERISH = new BiopsyTreatmentData("mmm", null, null, Lists.newArrayList());
    private final static BiopsyTreatmentData TREATMENT_WRONG_DRUG_DATA =
            new BiopsyTreatmentData("Yes", null, null, Lists.newArrayList(DRUG_NULL, DRUG_WRONG));

    private final static BiopsyTreatmentData TREATMENT_JAN_JAN =
            new BiopsyTreatmentData("Yes", JAN2015, JAN2015, Lists.newArrayList(DRUG_JAN_JAN));
    private final static BiopsyTreatmentData TREATMENT_JAN_FEB =
            new BiopsyTreatmentData("Yes", JAN2015, FEB2015, Lists.newArrayList(DRUG_JAN_FEB));
    private final static BiopsyTreatmentData TREATMENT_JAN_MAR =
            new BiopsyTreatmentData("Yes", JAN2015, MAR2015, Lists.newArrayList(DRUG_JAN_MAR));
    private final static BiopsyTreatmentData TREATMENT_JAN_ONGOING = new BiopsyTreatmentData("Yes", JAN2015, null, Lists.newArrayList());
    private final static BiopsyTreatmentData TREATMENT_FEB_ONGOING = new BiopsyTreatmentData("Yes", FEB2015, null, Lists.newArrayList());

    @Test
    public void reportsMissingDrugData() {
        final List<ValidationFinding> findings = PatientValidator.validateDrugData(CPCT_ID, DRUG_NULL);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_DRUG_START));
        assertTrue(findingsFields.contains(FIELD_DRUG));
    }

    @Test
    public void reportsIncorrectDrugData() {
        final List<ValidationFinding> findings = PatientValidator.validateDrugData(CPCT_ID, DRUG_WRONG);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(fields(FIELD_DRUG_START, FIELD_DRUG_END)));
        assertTrue(findingsFields.contains(FIELD_DRUG));
    }

    @Test
    public void reportsMissingTreatmentGiven() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentData(CPCT_ID, TREATMENT_GIVEN_NULL);
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_TREATMENT_GIVEN));
    }

    @Test
    public void reportsMissingTreatmentData() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentData(CPCT_ID, TREATMENT_GIVEN_EMPTY);
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FORM_TREATMENT));
    }

    @Test
    public void reportsTreatmentGivenNoWithData() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentData(CPCT_ID, TREATMENT_NOTGIVEN_DATA);
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_TREATMENT_GIVEN));
    }

    @Test
    public void reportsTreatmentGivenGibberish() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentData(CPCT_ID, TREATMENT_GIVEN_GIBBERISH);
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_TREATMENT_GIVEN));
    }

    @Test
    public void reportsDrugFindingsForTreatment() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentData(CPCT_ID, TREATMENT_WRONG_DRUG_DATA);
        assertEquals(4, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_DRUG_START));
        assertTrue(findingsFields.contains(FIELD_DRUG));
    }

    @Test
    public void reportsWrongTreatmentTimeline() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(CPCT_ID, Lists.newArrayList(TREATMENT_JAN_MAR, TREATMENT_JAN_JAN, TREATMENT_JAN_FEB));
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FORM_TREATMENT));
    }

    @Test
    public void reportsWrongTreatmentTimelineOngoing() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(CPCT_ID, Lists.newArrayList(TREATMENT_JAN_JAN, TREATMENT_JAN_ONGOING));
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FORM_TREATMENT));
    }

    @Test
    public void reportsTwoOngoingTreatments() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatments(CPCT_ID, Lists.newArrayList(TREATMENT_JAN_ONGOING, TREATMENT_FEB_ONGOING));
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FORM_TREATMENT));
    }

    @Test
    public void doesNotReportCorrectDeathTimeline() {
        final List<ValidationFinding> findings =
                PatientValidator.validateDeathDate(CPCT_ID, MAR2015, Lists.newArrayList(TREATMENT_JAN_JAN, TREATMENT_JAN_FEB));
        assertEquals(0, findings.size());
    }

    @Test
    public void reportsDeathDateBeforeEndOfTreatment() {
        final List<ValidationFinding> findings =
                PatientValidator.validateDeathDate(CPCT_ID, MAR2015, Lists.newArrayList(TREATMENT_JAN_ONGOING, TREATMENT_JAN_FEB));
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(fields(FIELD_DEATH_DATE, FORM_TREATMENT)));
    }
}
