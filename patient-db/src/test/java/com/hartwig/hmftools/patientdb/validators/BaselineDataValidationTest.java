package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.baselineBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocation;
import com.hartwig.hmftools.patientdb.data.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;

import org.junit.Test;

public class BaselineDataValidationTest {
    private static final String PATIENT_IDENTIFIER = "CPCT01020000";
    private static final String HOSPITAL = "Test Hospital";

    private static final BaselineData EMPTY_BASELINE = baselineBuilder().build();

    private static final BaselineData BASELINE_DATA_MISSING_LOCATION_MAPPING = baselineBuilder().hospital(HOSPITAL)
            .curatedTumorLocation(ImmutableCuratedTumorLocation.of(null, null, "some_location"))
            .build();

    @Test
    public void reportsMissingFields() {
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
}
