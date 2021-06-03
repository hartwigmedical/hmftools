package com.hartwig.hmftools.virusinterpreter.algo;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusBlacklistFileTest {

    private static final String VIRUS_BLACKLIST_TSV = Resources.getResource("viral_reporting/virus_blacklist.tsv").getPath();

    @Test
    public void canReadVirusBlacklistTsv() throws IOException {
        VirusBlacklistModel virusBlacklistModel = VirusBlacklistFile.buildFromTsv(VIRUS_BLACKLIST_TSV);
        assertEquals(1, virusBlacklistModel.genusCount());
        assertEquals(2, virusBlacklistModel.speciesCount());
    }
}