package com.hartwig.hmftools.common.clinical;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PatientPrimaryTumorFunctionsTest {

    @Test
    public void canFindPatientPrimaryTumors() {
        PatientPrimaryTumor cpct = createTestBuilder().patientIdentifier("CPCT02020202").location("loc1").type("type1").build();
        PatientPrimaryTumor colo = createTestBuilder().patientIdentifier("COLO829").location("loc2").type("type2").build();

        List<PatientPrimaryTumor> patientPrimaryTumors = Lists.newArrayList(cpct, colo);

        PatientPrimaryTumor primary1 = PatientPrimaryTumorFunctions.findPrimaryTumorForPatient(patientPrimaryTumors, "CPCT02020202");
        assertNotNull(primary1);
        assertEquals(cpct.location(), primary1.location());
        assertEquals(cpct.type(), primary1.type());

        PatientPrimaryTumor primary2 = PatientPrimaryTumorFunctions.findPrimaryTumorForPatient(patientPrimaryTumors, "COLO829");
        assertNotNull(primary2);
        assertEquals(colo.location(), primary2.location());
        assertEquals(colo.type(), primary2.type());

        PatientPrimaryTumor primary3 = PatientPrimaryTumorFunctions.findPrimaryTumorForPatient(patientPrimaryTumors, "any");
        assertNull(primary3);
    }

    @NotNull
    private static ImmutablePatientPrimaryTumor.Builder createTestBuilder() {
        return ImmutablePatientPrimaryTumor.builder()
                .subLocation("subLoc")
                .subType("subType")
                .extraDetails("details")
                .isOverridden(false);
    }
}