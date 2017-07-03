package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_BIRTH_YEAR1;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_BIRTH_YEAR2;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_BIRTH_YEAR3;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_ETHNICITY;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_PRIMARY_TUMOR_LOCATION;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_PRIMARY_TUMOR_LOCATION_OTHER;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_REGISTRATION_DATE;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_SEX;
import static com.hartwig.hmftools.patientdb.validators.PatientValidator.fields;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.data.PatientData;

import org.junit.Test;

public class PatientDataValidationTest {
    private final String CPCT_ID = "CPCT01020000";
    private final String HOSPITAL = "Test Hospital";
    private final PatientData PATIENT_DATA = new PatientData(CPCT_ID, null, null, null, HOSPITAL, null, null, null);

    @Test
    public void reportsMissingFields() {
        final List<ValidationFinding> findings = PatientValidator.validatePatientData(PATIENT_DATA);
        assertEquals(5, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_REGISTRATION_DATE));
        assertTrue(findingsFields.contains(FIELD_SEX));
        assertTrue(findingsFields.contains(FIELD_ETHNICITY));
        assertTrue(findingsFields.contains(fields(FIELD_BIRTH_YEAR1, FIELD_BIRTH_YEAR2, FIELD_BIRTH_YEAR3)));
        assertTrue(findingsFields.contains(fields(FIELD_PRIMARY_TUMOR_LOCATION, FIELD_PRIMARY_TUMOR_LOCATION_OTHER)));
    }
}
