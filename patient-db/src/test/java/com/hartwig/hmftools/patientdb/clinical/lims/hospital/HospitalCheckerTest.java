package com.hartwig.hmftools.patientdb.clinical.lims.hospital;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class HospitalCheckerTest {

    private static final String LIMS_DIRECTORY = Resources.getResource("lims").getPath();

    @Test
    public void canCreateFromLimsDirectoryWithoutWarnings() throws IOException {
        HospitalModel hospitalModel = HospitalModelFactory.fromLimsDirectory(LIMS_DIRECTORY);
        assertNotNull(hospitalModel);

        assertTrue(HospitalChecker.validateModelIntegrity(hospitalModel));
    }

    @Test
    public void canValidateIntegrityCorrectly() {
        HospitalModel emptyModel = ImmutableHospitalModel.builder().build();
        assertTrue(HospitalChecker.validateModelIntegrity(emptyModel));

        HospitalModel noAddressModel = ImmutableHospitalModel.builder()
                .putHospitalPersonsCPCT("HOSP",
                        ImmutableHospitalPersons.builder()
                                .hospitalPI("Test")
                                .requesterName("Test")
                                .requesterEmail("Test")
                                .build())
                .build();

        assertFalse(HospitalChecker.validateModelIntegrity(noAddressModel));
    }
}