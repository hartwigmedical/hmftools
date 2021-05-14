package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusBlacklistFileTest {

    private static final String VIRUS_BLACKLIST_TSV = Resources.getResource("virusbreakend/virus_blacklist.tsv").getPath();

    @Test
    public void readVirusBlacklistTsv() throws IOException {
        VirusBlackListModel2 virusBlacklistModel = VirusBlacklistFile.buildFromTsv(VIRUS_BLACKLIST_TSV);
        assertEquals(3, virusBlacklistModel.virusBlacklistCount());

        assertEquals("taxid_genus", virusBlacklistModel.checkTaxusForId(1));
        assertEquals("taxid_species", virusBlacklistModel.checkTaxusForId(2));

        assertTrue(virusBlacklistModel.checkVirusForBlacklisting(1));
        assertTrue(virusBlacklistModel.checkVirusForBlacklisting(2));
        assertFalse(virusBlacklistModel.checkVirusForBlacklisting(40));
    }
}