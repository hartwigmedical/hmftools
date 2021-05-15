package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusSummaryFileTest {

    private static final String VIRUS_SUMMARY_TSV = Resources.getResource("virusbreakend/virus_summary.tsv").getPath();

    @Test
    public void readVirusDbTsv() throws IOException {
        VirusSummaryModel virusSummaryModel = VirusSummaryFile.buildFromTsv(VIRUS_SUMMARY_TSV);
        assertEquals(1, virusSummaryModel.virusCount());

        assertTrue(virusSummaryModel.mapIdToVirusName(1));
        assertFalse(virusSummaryModel.mapIdToVirusName(2));

        assertEquals("virus", virusSummaryModel.findVirusSummary(1));
        assertNotEquals("virus1", virusSummaryModel.findVirusSummary(2));
    }
}