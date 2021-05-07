package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.*;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusDbFileTest {

    private static final String VIRUSDB_TSV = Resources.getResource("virusbreakend/virusdb.tsv").getPath();

    @Test
    public void readVirusDbTsv() throws IOException {
        VirusDbModel virusDbModel = VirusDbFile.buildFromTsv(VIRUSDB_TSV);
        assertEquals(1, virusDbModel.viruscount());

        assertTrue(virusDbModel.mapIdtoVirusName(1));
        assertFalse(virusDbModel.mapIdtoVirusName(2));
    }
}