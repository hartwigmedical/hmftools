package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusDbFileTest {

    private static final String VIRUSDB_TSV = Resources.getResource("virusbreakend/virusdb.tsv").getPath();

    @Test
    public void readVirusDbTsv() throws IOException {
        VirusDbModel virusDbModel = VirusDbFile.buildFromTsv(VIRUSDB_TSV);
        assertEquals(1, virusDbModel.virusCount());

        assertTrue(virusDbModel.mapIdToVirusName(1));
        assertFalse(virusDbModel.mapIdToVirusName(2));
    }
}