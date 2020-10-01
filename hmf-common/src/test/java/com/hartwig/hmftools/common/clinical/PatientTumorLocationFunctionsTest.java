package com.hartwig.hmftools.common.clinical;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class PatientTumorLocationFunctionsTest {

    @Test
    public void canFindPatientTumorLocations() {
        PatientTumorLocation cpct = ImmutablePatientTumorLocation.builder()
                .patientIdentifier("CPCT02020202")
                .primaryTumorLocation("tum1")
                .cancerSubtype("sub1")
                .build();

        PatientTumorLocation colo = ImmutablePatientTumorLocation.builder()
                .patientIdentifier("COLO829")
                .primaryTumorLocation("tum2")
                .cancerSubtype("sub2")
                .build();

        List<PatientTumorLocation> patientTumorLocations = Lists.newArrayList(cpct, colo);

        PatientTumorLocation location1 = PatientTumorLocationFunctions.findTumorLocationForSample(patientTumorLocations, "CPCT02020202T");
        assertNotNull(location1);
        assertEquals(cpct.primaryTumorLocation(), location1.primaryTumorLocation());
        assertEquals(cpct.cancerSubtype(), location1.cancerSubtype());

        PatientTumorLocation location2 = PatientTumorLocationFunctions.findTumorLocationForSample(patientTumorLocations, "COLO829T");
        assertNotNull(location2);
        assertEquals(colo.primaryTumorLocation(), location2.primaryTumorLocation());
        assertEquals(colo.cancerSubtype(), location2.cancerSubtype());

        PatientTumorLocation location3 = PatientTumorLocationFunctions.findTumorLocationForSample(patientTumorLocations, "any");
        assertNull(location3);
    }
}