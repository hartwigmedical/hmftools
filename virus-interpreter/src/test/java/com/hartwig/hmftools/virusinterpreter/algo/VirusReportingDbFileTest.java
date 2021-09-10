package com.hartwig.hmftools.virusinterpreter.algo;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusReportingDbFileTest {

    private static final String VIRUS_REPORTING_DB_TSV = Resources.getResource("virus_interpreter/virus_reporting_db.tsv").getPath();

    @Test
    public void canReadVirusReportingDbTsv() throws IOException {
        VirusReportingDbModel virusReportingDbModel = VirusReportingDbFile.buildFromTsv(VIRUS_REPORTING_DB_TSV);
        assertEquals(1, virusReportingDbModel.count());

        assertTrue(virusReportingDbModel.hasInterpretation(1));
        assertFalse(virusReportingDbModel.hasInterpretation(2));

        assertEquals("MCV", virusReportingDbModel.interpretVirusSpecies(1));
        assertEquals(Integer.valueOf(90), virusReportingDbModel.nonIntegratedMinimalCoverage(1));
        assertNull(virusReportingDbModel.integratedMinimalCoverage(1));

        assertNotEquals("HPV", virusReportingDbModel.interpretVirusSpecies(2));
    }
}