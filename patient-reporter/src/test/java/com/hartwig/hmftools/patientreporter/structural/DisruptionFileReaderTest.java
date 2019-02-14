package com.hartwig.hmftools.patientreporter.structural;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class DisruptionFileReaderTest {

    private static final String DISRUPTION_FILE = Resources.getResource("test_run/svAnalysis/CPCT11111111T_disruptions.csv").getPath();

    @Test
    public void canReadFromFile() throws IOException {
        List<Disruption> disruptions = DisruptionFileReader.fromDisruptionFile(DISRUPTION_FILE);

        assertEquals(3, disruptions.size());
        assertNotNull(disruptions.get(0).ploidy());
        // 2nd disruption has missing ploidy
        assertNull(disruptions.get(1).ploidy());
        // 3rd disruption has ploidy 0.00000
        assertNull(disruptions.get(2).ploidy());
    }
}
