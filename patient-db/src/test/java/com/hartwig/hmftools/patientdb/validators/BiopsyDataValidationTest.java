package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.baselineBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyBuilder;

import static org.junit.Assert.assertEquals;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.data.BiopsyData;

import org.junit.Test;

public class BiopsyDataValidationTest {
    private static final String PATIENT_IDENTIFIER = "CPCT01020000";

    private static final LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private static final LocalDate MAR2016 = LocalDate.parse("2016-03-01");

    private static final BiopsyData BIOPSY_NULL = biopsyBuilder().date(null).build();
    private static final BiopsyData BIOPSY_FEB1 = biopsyBuilder().date(FEB2015).site("1").location("").build();
    private static final BiopsyData BIOPSY_FEB2 = biopsyBuilder().date(FEB2015).site("2").location("").build();

    @Test
    public void reportsMissingFields() {
        final List<ValidationFinding> findings = PatientValidator.validateBiopsyData(PATIENT_IDENTIFIER, BIOPSY_NULL);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void reportsAllFieldsEmpty() {
        final List<ValidationFinding> findings =
                PatientValidator.validateBiopsies(PATIENT_IDENTIFIER, Lists.newArrayList(BIOPSY_NULL, BIOPSY_FEB1, BIOPSY_FEB2));
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
    }

    @Test
    public void reportsBiopsyBeforeInformedConsent() {
        final List<ValidationFinding> findings = PatientValidator.validateInformedConsentDate(PATIENT_IDENTIFIER,
                baselineBuilder().informedConsentDate(MAR2016).build(),
                Lists.newArrayList(BIOPSY_FEB1));
        assertEquals(1, findings.size());
    }
}
