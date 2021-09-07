package com.hartwig.hmftools.virusinterpreter.algo;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusWhitelistFileTest {

    private static final String VIRUS_WHITELIST_TSV = Resources.getResource("virus_interpreter/virus_whitelist.tsv").getPath();

    @Test
    public void canReadVirusWhitelistTsv() throws IOException {
        VirusWhitelistModel virusWhitelistModel = VirusWhitelistFile.buildFromTsv(VIRUS_WHITELIST_TSV);
        assertEquals(1, virusWhitelistModel.count());

        assertTrue(virusWhitelistModel.hasInterpretation(1));
        assertFalse(virusWhitelistModel.hasInterpretation(2));

        assertEquals("MCV", virusWhitelistModel.interpretVirusSpecies(1));
        assertTrue(virusWhitelistModel.displayVirusOnSummaryReport(1));
        assertEquals(Integer.valueOf(90), virusWhitelistModel.nonIntegratedMinimalCoverage(1));
        assertNull(virusWhitelistModel.integratedMinimalCoverage(1));

        assertNotEquals("HPV", virusWhitelistModel.interpretVirusSpecies(2));
    }
}