package com.hartwig.hmftools.common.hospital;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class CenterModelFactoryTest {

    private static final String CENTER_RESOURCE = Resources.getResource("hospital/centers.csv").getPath();
    private static final String CENTER_RESOURCE_MANUAL = Resources.getResource("hospital/manual_mapping.csv").getPath();

    @Test
    public void canReadFromFile() throws IOException {
        final CenterModel centerModel = CenterModelFactory.readFromCSV(CENTER_RESOURCE, CENTER_RESOURCE_MANUAL);
        assertEquals(2, centerModel.centerPerId().size());
        assertEquals(2, centerModel.centerPerHospital().size());
    }
}
