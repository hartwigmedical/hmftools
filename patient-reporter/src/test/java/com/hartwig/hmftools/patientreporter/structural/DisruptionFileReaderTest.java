package com.hartwig.hmftools.patientreporter.structural;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class DisruptionFileReaderTest {

    private static final String DISRUPTION_FILE = Resources.getResource("test_run/svAnalysis/CPCT11111111T_disruptions.csv").getPath();

    @Test
    public void canReadFromFile() throws IOException {
        assertEquals(1, DisruptionFileReader.fromDisruptionFile(DISRUPTION_FILE).size());
    }
}
