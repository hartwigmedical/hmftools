package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.baselineBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyTreatmentBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.drugBuilder;
import static com.hartwig.hmftools.patientdb.readers.cpct.BaselineReader.FIELD_INFORMED_CONSENT_DATE;
import static com.hartwig.hmftools.patientdb.readers.cpct.BiopsyReader.FIELD_BIOPSY_DATE;
import static com.hartwig.hmftools.patientdb.readers.cpct.BiopsyReader.FIELD_LOCATION;
import static com.hartwig.hmftools.patientdb.readers.cpct.BiopsyReader.FIELD_SITE;
import static com.hartwig.hmftools.patientdb.readers.cpct.BiopsyReader.FIELD_SITE_OTHER;
import static com.hartwig.hmftools.patientdb.readers.cpct.BiopsyReader.FORM_BIOPS;
import static com.hartwig.hmftools.patientdb.readers.cpct.BiopsyTreatmentReader.FORM_TREATMENT;
import static com.hartwig.hmftools.patientdb.validators.PatientValidator.fields;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;

import org.junit.Test;

public class BiopsyDataValidationTest {
    private static final String PATIENT_IDENTIFIER = "CPCT01020000";
    private static final LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private static final LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private static final LocalDate MAR2016 = LocalDate.parse("2016-03-01");

    private static final BiopsyData BIOPSY_NULL = biopsyBuilder().date(null).build();
    private static final BiopsyData BIOPSY_FEB1 = biopsyBuilder().date(FEB2015).site("1").location("").build();

    private static final BiopsyData BIOPSY_FEB2 = biopsyBuilder().date(FEB2015).site("2").location("").build();
    private static final BiopsyTreatmentData TREATMENT_JAN_FEB =
            biopsyTreatmentBuilder().addDrugs(drugBuilder().startDate(JAN2015).endDate(FEB2015).build()).build();

    @Test
    public void reportsMissingFields() {
        final List<ValidationFinding> findings = PatientValidator.validateBiopsyData(PATIENT_IDENTIFIER, BIOPSY_NULL);
        assertEquals(3, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_BIOPSY_DATE));
        assertTrue(findingsFields.contains(fields(FIELD_SITE, FIELD_SITE_OTHER)));
        assertTrue(findingsFields.contains(fields(FIELD_LOCATION)));
    }

    @Test
    public void reportsBiopsiesEmpty() {
        final List<ValidationFinding> findings =
                PatientValidator.validateBiopsies(PATIENT_IDENTIFIER, Lists.newArrayList(), Lists.newArrayList());
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FORM_BIOPS));
    }

    @Test
    public void reportsAllFieldsEmpty() {
        final List<ValidationFinding> findings = PatientValidator.validateBiopsies(PATIENT_IDENTIFIER,
                Lists.newArrayList(BIOPSY_NULL, BIOPSY_FEB1, BIOPSY_FEB2),
                Lists.newArrayList());
        assertEquals(3, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_BIOPSY_DATE));
        assertTrue(findingsFields.contains(fields(FIELD_SITE, FIELD_SITE_OTHER)));
        assertTrue(findingsFields.contains(fields(FIELD_LOCATION)));
    }

    @Test
    public void reportsBiopsyBeforeInformedConsent() {
        final List<ValidationFinding> findings = PatientValidator.validateInformedConsentDate(PATIENT_IDENTIFIER,
                baselineBuilder().informedConsentDate(MAR2016).build(),
                Lists.newArrayList(BIOPSY_FEB1));
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(fields(FIELD_INFORMED_CONSENT_DATE, FIELD_BIOPSY_DATE)));
    }

    @Test
    public void reportsTreatmentBeforeBiopsy() {
        final List<ValidationFinding> findings = PatientValidator.validateBiopsies(PATIENT_IDENTIFIER,
                Lists.newArrayList(BIOPSY_FEB1),
                Lists.newArrayList(TREATMENT_JAN_FEB));
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(fields(FORM_TREATMENT, FORM_BIOPS)));
    }
}
