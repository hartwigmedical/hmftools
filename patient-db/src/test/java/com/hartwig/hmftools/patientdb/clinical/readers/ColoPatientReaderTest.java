package com.hartwig.hmftools.patientdb.clinical.readers;

import static org.junit.Assert.assertNotNull;

import org.junit.Test;

public class ColoPatientReaderTest {

    @Test
    public void canReadColo() {
        assertNotNull(new ColoPatientReader().read("COLO829T"));
    }
}