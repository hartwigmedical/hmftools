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
        PatientPrimaryTumor cpct =
                createTestBuilder().patientIdentifier("CPCT02020202").primaryTumorLocation("loc1").primaryTumorType("type1").build();

        PatientPrimaryTumor colo =
                createTestBuilder().patientIdentifier("COLO829").primaryTumorLocation("loc2").primaryTumorType("type2").build();

        List<PatientPrimaryTumor> patientPrimaryTumors = Lists.newArrayList(cpct, colo);

        PatientPrimaryTumor primary1 =
                PatientPrimaryTumorFunctions.findPrimaryTumorForSample(patientPrimaryTumors, "CPCT02020202T");
        assertNotNull(primary1);
        assertEquals(cpct.primaryTumorLocation(), primary1.primaryTumorLocation());
        assertEquals(cpct.primaryTumorType(), primary1.primaryTumorType());

        PatientPrimaryTumor primary2 = PatientPrimaryTumorFunctions.findPrimaryTumorForSample(patientPrimaryTumors, "COLO829T");
        assertNotNull(primary2);
        assertEquals(colo.primaryTumorLocation(), primary2.primaryTumorLocation());
        assertEquals(colo.primaryTumorType(), primary2.primaryTumorType());

        PatientPrimaryTumor primary3 = PatientPrimaryTumorFunctions.findPrimaryTumorForSample(patientPrimaryTumors, "any");
        assertNull(primary3);
    }

    @NotNull
    private static ImmutablePatientPrimaryTumor.Builder createTestBuilder() {
        return ImmutablePatientPrimaryTumor.builder()
                .primaryTumorSubLocation("subLoc")
                .primaryTumorSubType("subType")
                .primaryTumorExtraDetails("details")
                .doids(Lists.newArrayList("doid"))
                .isOverridden(false);
    }
}