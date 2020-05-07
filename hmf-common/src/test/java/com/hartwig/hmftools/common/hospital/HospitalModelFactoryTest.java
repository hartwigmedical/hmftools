package com.hartwig.hmftools.common.hospital;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class HospitalModelFactoryTest {

    private static final String HOSPITAL_RESOURCE = Resources.getResource("hospital").getPath();

    @Test
    public void canReadFromHospitalDirectory() throws IOException {
        HospitalModel hospitalModel = HospitalModelFactory.fromHospitalDirectory(HOSPITAL_RESOURCE);

       // assertEquals(2, hospitalModel.hospitalCount());
        assertEquals(1, hospitalModel.sampleHospitalMapping().size());
     //   assertEquals(1, hospitalModel.hospitalCoreMap().size());
    }
}
