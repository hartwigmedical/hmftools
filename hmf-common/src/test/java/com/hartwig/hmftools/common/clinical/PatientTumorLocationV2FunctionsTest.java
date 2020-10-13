package com.hartwig.hmftools.common.clinical;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PatientTumorLocationV2FunctionsTest {

    @Test
    public void canFindPatientTumorLocations() {
        PatientTumorLocationV2 cpct =
                createTestBuilder().patientIdentifier("CPCT02020202").primaryTumorLocation("loc1").primaryTumorType("type1").build();

        PatientTumorLocationV2 colo =
                createTestBuilder().patientIdentifier("COLO829").primaryTumorLocation("loc2").primaryTumorType("type2").build();

        List<PatientTumorLocationV2> patientTumorLocations = Lists.newArrayList(cpct, colo);

        PatientTumorLocationV2 location1 =
                PatientTumorLocationV2Functions.findTumorLocationForSample(patientTumorLocations, "CPCT02020202T");
        assertNotNull(location1);
        assertEquals(cpct.primaryTumorLocation(), location1.primaryTumorLocation());
        assertEquals(cpct.primaryTumorType(), location1.primaryTumorType());

        PatientTumorLocationV2 location2 = PatientTumorLocationV2Functions.findTumorLocationForSample(patientTumorLocations, "COLO829T");
        assertNotNull(location2);
        assertEquals(colo.primaryTumorLocation(), location2.primaryTumorLocation());
        assertEquals(colo.primaryTumorType(), location2.primaryTumorType());

        PatientTumorLocationV2 location3 = PatientTumorLocationV2Functions.findTumorLocationForSample(patientTumorLocations, "any");
        assertNull(location3);
    }

    @NotNull
    private static ImmutablePatientTumorLocationV2.Builder createTestBuilder() {
        return ImmutablePatientTumorLocationV2.builder()
                .primaryTumorSubLocation("subLoc")
                .primaryTumorSubType("subType")
                .primaryTumorExtraDetails("details")
                .doids(Lists.newArrayList("doid"))
                .isOverridden(false);
    }
}