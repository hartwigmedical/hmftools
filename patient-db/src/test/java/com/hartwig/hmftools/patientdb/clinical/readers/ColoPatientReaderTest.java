package com.hartwig.hmftools.patientdb.clinical.readers;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.patientdb.clinical.curators.CuratorTestFactory;

import org.junit.Test;

public class ColoPatientReaderTest {

    @Test
    public void canReadColo() {
        assertNotNull(new ColoPatientReader(CuratorTestFactory.primaryTumorCurator()).read("COLO829T", "Melanoma"));
    }
}