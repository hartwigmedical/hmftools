package com.hartwig.hmftools.virusinterpreter.algo;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusBlacklistingDbFileTest
{
    private static final String VIRUS_BLACKLISTING_DB_TSV = Resources.getResource("virus_interpreter/virus_blacklisting_db.tsv").getPath();

    @Test
    public void canReadVirusBlacklistingDbTsv() throws IOException
    {
        List<Integer> blacklistedTaxids = VirusBlacklistingDbFile.loadFromTsv(VIRUS_BLACKLISTING_DB_TSV);
        assertEquals(3, blacklistedTaxids.size());
    }
}
