package com.hartwig.hmftools.patientdb.clinical.readers;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.patientdb.clinical.curators.TestCuratorFactory;

import org.junit.Test;

public class ColoPatientReaderTest {

    @Test
    public void canReadColo() {
        assertNotNull(new ColoPatientReader(TestCuratorFactory.primaryTumorCurator()).read("COLO829T", "Melanoma"));
    }
}