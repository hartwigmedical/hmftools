package com.hartwig.hmftools.common.clinical;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PatientTumorLocationFunctionsTest {

    @Test
    public void canFindPatientTumorLocations() {
        PatientTumorLocation cpct =
                createTestBuilder().patientIdentifier("CPCT02020202").primaryTumorLocation("loc1").primaryTumorType("type1").build();

        PatientTumorLocation colo =
                createTestBuilder().patientIdentifier("COLO829").primaryTumorLocation("loc2").primaryTumorType("type2").build();

        List<PatientTumorLocation> patientTumorLocations = Lists.newArrayList(cpct, colo);

        PatientTumorLocation location1 =
                PatientTumorLocationFunctions.findTumorLocationForSample(patientTumorLocations, "CPCT02020202T");
        assertNotNull(location1);
        assertEquals(cpct.primaryTumorLocation(), location1.primaryTumorLocation());
        assertEquals(cpct.primaryTumorType(), location1.primaryTumorType());

        PatientTumorLocation location2 = PatientTumorLocationFunctions.findTumorLocationForSample(patientTumorLocations, "COLO829T");
        assertNotNull(location2);
        assertEquals(colo.primaryTumorLocation(), location2.primaryTumorLocation());
        assertEquals(colo.primaryTumorType(), location2.primaryTumorType());

        PatientTumorLocation location3 = PatientTumorLocationFunctions.findTumorLocationForSample(patientTumorLocations, "any");
        assertNull(location3);
    }

    @NotNull
    private static ImmutablePatientTumorLocation.Builder createTestBuilder() {
        return ImmutablePatientTumorLocation.builder()
                .primaryTumorSubLocation("subLoc")
                .primaryTumorSubType("subType")
                .primaryTumorExtraDetails("details")
                .doids(Lists.newArrayList("doid"))
                .isOverridden(false);
    }
}