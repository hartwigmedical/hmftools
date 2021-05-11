package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.*;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusBlacklistFileTest {

    private static final String VIRUS_BLACKLIST_TSV = Resources.getResource("virusbreakend/virus_blacklist.tsv").getPath();

    @Test
    public void readVirusBlacklistTsv() throws IOException {
        VirusBlackListModel virusBlacklistModel = VirusBlacklistFile.buildFromTsv(VIRUS_BLACKLIST_TSV);
        assertEquals(3, virusBlacklistModel.virusBlacklistcount());

        assertTrue(virusBlacklistModel.checkTaxus("taxid_genus"));
        assertTrue(virusBlacklistModel.checkTaxus("taxid_species"));
        assertFalse(virusBlacklistModel.checkTaxus("species"));

        assertTrue(virusBlacklistModel.checkVirusForBlacklisting(1));
        assertTrue(virusBlacklistModel.checkVirusForBlacklisting(2));
        assertFalse(virusBlacklistModel.checkVirusForBlacklisting(40));

    }
}