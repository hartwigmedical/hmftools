package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_BIRTH_YEAR1;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_BIRTH_YEAR2;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_BIRTH_YEAR3;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_HOSPITAL1;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_HOSPITAL2;
import static com.hartwig.hmftools.patientdb.readers.CpctPatientReader.FIELD_INFORMED_CONSENT_DATE;
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
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutableBaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedCancerType;
import com.hartwig.hmftools.patientdb.data.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;

import org.junit.Test;

public class BaselineDataValidationTest {
    private static final String PATIENT_IDENTIFIER = "CPCT01020000";
    private static final String HOSPITAL = "Test Hospital";

    private static final BaselineData EMPTY_BASELINE = ImmutableBaselineData.builder().build();

    private static final BaselineData BASELINE_DATA_MISSING_LOCATION_MAPPING = ImmutableBaselineData.builder()
            .hospital(HOSPITAL)
            .cancerType(ImmutableCuratedCancerType.of(null, null, "some_location"))
            .build();

    @Test
    public void reportsMissingFields() {
        final List<ValidationFinding> findings = PatientValidator.validateBaselineData(PATIENT_IDENTIFIER, EMPTY_BASELINE);
        assertEquals(6, findings.size());
        findings.stream().map(ValidationFinding::patientId).forEach(id -> assertEquals(PATIENT_IDENTIFIER, id));
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());

        assertTrue(findingsFields.contains(fields(FIELD_REGISTRATION_DATE1, FIELD_REGISTRATION_DATE2)));
        assertTrue(findingsFields.contains(FIELD_INFORMED_CONSENT_DATE));
        assertTrue(findingsFields.contains(FIELD_SEX));
        assertTrue(findingsFields.contains(fields(FIELD_BIRTH_YEAR1, FIELD_BIRTH_YEAR2, FIELD_BIRTH_YEAR3)));
        assertTrue(findingsFields.contains(fields(FIELD_PRIMARY_TUMOR_LOCATION, FIELD_PRIMARY_TUMOR_LOCATION_OTHER)));
        assertTrue(findingsFields.contains(fields(FIELD_HOSPITAL1, FIELD_HOSPITAL2)));
    }

    @Test
    public void reportMissingPreTreatment() {
        PreTreatmentData emptyData = ImmutablePreTreatmentData.builder().formLocked(false).formStatus(FormStatusState.UNKNOWN).build();

        assertEquals(2, PatientValidator.validatePreTreatmentData(PATIENT_IDENTIFIER, emptyData).size());

        PreTreatmentData actualData = ImmutablePreTreatmentData.builder()
                .radiotherapyGiven("Yes")
                .treatmentGiven("No")
                .formLocked(false)
                .formStatus(FormStatusState.UNKNOWN)
                .build();
        assertEquals(0, PatientValidator.validatePreTreatmentData(PATIENT_IDENTIFIER, actualData).size());

        PreTreatmentData onlyRadioTherapyPresent =
                ImmutablePreTreatmentData.builder().radiotherapyGiven("Yes").formLocked(false).formStatus(FormStatusState.UNKNOWN).build();
        assertEquals(1, PatientValidator.validatePreTreatmentData(PATIENT_IDENTIFIER, onlyRadioTherapyPresent).size());

        PreTreatmentData onlyTreatmentPresent =
                ImmutablePreTreatmentData.builder().treatmentGiven("Yes").formLocked(false).formStatus(FormStatusState.UNKNOWN).build();
        assertEquals(1, PatientValidator.validatePreTreatmentData(PATIENT_IDENTIFIER, onlyTreatmentPresent).size());
    }

    @Test
    public void reportsFailureToCuratePrimaryTumorLocation() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTumorLocationCuration(PATIENT_IDENTIFIER, BASELINE_DATA_MISSING_LOCATION_MAPPING);
        assertEquals(1, findings.size());
        assertEquals("tumorLocationCuration", findings.get(0).level());
        final List<String> findingsFields = findings.stream().map(ValidationFinding::ecrfItem).collect(Collectors.toList());
        assertTrue(findingsFields.contains(fields(FIELD_PRIMARY_TUMOR_LOCATION, FIELD_PRIMARY_TUMOR_LOCATION_OTHER)));
    }
}
