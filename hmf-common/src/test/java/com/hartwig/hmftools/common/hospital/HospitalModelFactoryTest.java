package com.hartwig.hmftools.common.hospital;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class HospitalModelFactoryTest {

    private static final String HOSPITAL_RESOURCE = Resources.getResource("hospital/hospitals.csv").getPath();
    private static final String HOSPITAL_SAMPLE_MAPPING_RESOURCE = Resources.getResource("hospital/sample_hospital_mapping.csv").getPath();

    @Test
    public void canReadFromFile() throws IOException {
        final HospitalModel hospitalModel = HospitalModelFactory.readFromCSV(HOSPITAL_RESOURCE, HOSPITAL_SAMPLE_MAPPING_RESOURCE);
        assertEquals(2, hospitalModel.hospitalPerId().size());
        assertEquals(2, hospitalModel.hospitalPerHospital().size());
    }
}
