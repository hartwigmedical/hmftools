package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_BIRTH_YEAR1;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_BIRTH_YEAR2;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_BIRTH_YEAR3;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_PRIMARY_TUMOR_LOCATION;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_PRIMARY_TUMOR_LOCATION_OTHER;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_REGISTRATION_DATE1;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_REGISTRATION_DATE2;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_SEX;
import static com.hartwig.hmftools.patientdb.validators.PatientValidator.fields;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocation;
import com.hartwig.hmftools.patientdb.data.ImmutablePatientData;
import com.hartwig.hmftools.patientdb.data.PatientData;

import org.junit.Test;

public class PatientDataValidationTest {
    private final String CPCT_ID = "CPCT01020000";
    private final String HOSPITAL = "Test Hospital";
    private final PatientData PATIENT_DATA =
            ImmutablePatientData.of(CPCT_ID, null, null, HOSPITAL, null, ImmutableCuratedTumorLocation.of(null, null, null), null,
                    FormStatusState.UNKNOWN, false, FormStatusState.UNKNOWN, false, FormStatusState.UNKNOWN, false, FormStatusState.UNKNOWN,
                    false, FormStatusState.UNKNOWN, false);
    private final PatientData PATIENT_DATA_MISSING_LOCATION_MAPPING =
            ImmutablePatientData.of(CPCT_ID, null, null, HOSPITAL, null, ImmutableCuratedTumorLocation.of(null, null, "some_location"),
                    null, FormStatusState.UNKNOWN, false, FormStatusState.UNKNOWN, false, FormStatusState.UNKNOWN, false,
                    FormStatusState.UNKNOWN, false, FormStatusState.UNKNOWN, false);

    @Test
    public void reportsMissingFields() {
        final List<ValidationFinding> findings = PatientValidator.validatePatientData(PATIENT_DATA);
        assertEquals(6, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(CPCT_ID, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(FIELD_REGISTRATION_DATE2));
        assertTrue(findingsFields.contains(FIELD_REGISTRATION_DATE1));
        assertTrue(findingsFields.contains(FIELD_SEX));
        assertTrue(findingsFields.contains(FIELD_BIRTH_YEAR1));
        assertTrue(findingsFields.contains(fields(FIELD_BIRTH_YEAR2, FIELD_BIRTH_YEAR3)));
        assertTrue(findingsFields.contains(fields(FIELD_PRIMARY_TUMOR_LOCATION, FIELD_PRIMARY_TUMOR_LOCATION_OTHER)));
    }

    @Test
    public void reportsFailureToCuratePrimaryTumorLocation() {
        final List<ValidationFinding> findings = PatientValidator.validateTumorLocationCuration(PATIENT_DATA_MISSING_LOCATION_MAPPING);
        assertEquals(1, findings.size());
        assertEquals("tumorLocationCuration", findings.get(0).level());
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(fields(FIELD_PRIMARY_TUMOR_LOCATION, FIELD_PRIMARY_TUMOR_LOCATION_OTHER)));
    }
}
