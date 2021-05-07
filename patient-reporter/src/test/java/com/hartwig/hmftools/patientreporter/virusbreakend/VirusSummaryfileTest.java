package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.*;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusSummaryfileTest {

    private static final String VIRUS_SUMMARY_TSV = Resources.getResource("virusbreakend/virus_summary.tsv").getPath();

    @Test
    public void readVirusDbTsv() throws IOException {
        VirusSummaryModel virusSummaryModel = VirusSummaryfile.buildFromTsv(VIRUS_SUMMARY_TSV);
        assertEquals(1, virusSummaryModel.viruscount());

        assertTrue(virusSummaryModel.mapIdtoVirusName(1));
        assertFalse(virusSummaryModel.mapIdtoVirusName(2));

        assertEquals("virus", virusSummaryModel.findVirusSummary(1));
        assertNotEquals("virus1", virusSummaryModel.findVirusSummary(2));

    }

}