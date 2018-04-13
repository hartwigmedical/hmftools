package com.hartwig.hmftools.patientdb.validators;

import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyTreatmentBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.biopsyTreatmentResponseBuilder;
import static com.hartwig.hmftools.patientdb.data.TestDatamodelFactory.drugBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.ImmutableBiopsyTreatmentResponseData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class TreatmentResponseValidationTest {
    private static final String CPCT_ID = "CPCT01020000";
    private static final LocalDate JAN2015 = LocalDate.parse("2015-01-01");
    private static final LocalDate FEB2015 = LocalDate.parse("2015-02-01");
    private static final LocalDate MAR2015 = LocalDate.parse("2015-03-01");

    private static final BiopsyTreatmentData TREATMENT_JAN_MAR =
            biopsyTreatmentBuilder().addDrugs(drugWithStartAndEndDate(JAN2015, MAR2015)).build();
    private static final BiopsyTreatmentData TREATMENT_JAN_ONGOING =
            biopsyTreatmentBuilder().addDrugs(drugWithStartAndEndDate(JAN2015, null)).build();
    private static final BiopsyTreatmentData TREATMENT_JAN_JAN =
            biopsyTreatmentBuilder().addDrugs(drugWithStartAndEndDate(JAN2015, JAN2015)).build();

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
    public void reportsMeasurementDoneNull() {
        final List<ValidationFinding> findings = validateNonFirstResponse(CPCT_ID, RESPONSE_NULL);
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsOnlyResponseFilledIn() {
        final List<ValidationFinding> findings = validateNonFirstResponse(CPCT_ID, RESPONSE_ONLY);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(CPCT_ID, id));
    }

    @Test
    public void reportsMeasurementDoneMissingData() {
        final List<ValidationFinding> findings = validateNonFirstResponse(CPCT_ID, RESPONSE_MISSING_DATA);
        assertEquals(2, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(CPCT_ID, id));
    }

    @Test
    public void acceptMissingDataForFirstResponse() {
        final List<ValidationFinding> findings = PatientValidator.validateTreatmentResponse(CPCT_ID, RESPONSE_MISSING_DATA, true);
        assertEquals(1, findings.size());
        findings.stream().map(ValidationFinding::patientIdentifier).forEach(id -> assertEquals(CPCT_ID, id));
    }

    @Test
    public void reportsMeasurementDoneNoWithData() {
        final List<ValidationFinding> findings = validateNonFirstResponse(CPCT_ID, RESPONSE_MEASUREMENT_NO_WITH_DATA);
        assertEquals(2, findings.size());
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
    }

    @Test
    public void reportsMissingTreatmentData() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatmentResponses(CPCT_ID, Lists.newArrayList(), Lists.newArrayList(RESPONSE_JAN2015));
        assertEquals(1, findings.size());
    }

    @Test
    public void reportsMissingResponseForTreatmentMoreThan16Weeks() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatmentResponses(CPCT_ID, Lists.newArrayList(TREATMENT_JAN_ONGOING), Lists.newArrayList());
        assertEquals(1, findings.size());
    }

    @Test
    public void doesNotReportMissingResponseForTreatmentEndedImmediately() {
        final List<ValidationFinding> findings =
                PatientValidator.validateTreatmentResponses(CPCT_ID, Lists.newArrayList(TREATMENT_JAN_JAN), Lists.newArrayList());
        assertEquals(0, findings.size());
    }

    @NotNull
    private static DrugData drugWithStartAndEndDate(@Nullable LocalDate startDate, @Nullable LocalDate endDate) {
        return drugBuilder().startDate(startDate).endDate(endDate).build();
    }

    @NotNull
    private static List<ValidationFinding> validateNonFirstResponse(@NotNull String patientId,
            @NotNull BiopsyTreatmentResponseData response) {
        return PatientValidator.validateTreatmentResponse(patientId, response, false);
    }
}
