package com.hartwig.hmftools.boggs.flagstatreader;

import com.google.common.io.Resources;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import static org.junit.Assert.assertEquals;

public class SambambaFlagstatParserTest {

    @Test
    public void canParseExampleFile() throws IOException {
        URL exampleFlagstatURL = Resources.getResource("flagstats/example.flagstat");
        String exampleFlagstatFile = exampleFlagstatURL.getPath();

        FlagstatParser parser = new SambambaFlagstatParser();
        Flagstat flagstat = parser.parse(new File(exampleFlagstatFile));

        assertEquals(0, flagstat.qcPassedReads().getTotal());
        assertEquals(1, flagstat.qcPassedReads().getSecondary());
        assertEquals(2, flagstat.qcPassedReads().getSupplementary());
        assertEquals(3, flagstat.qcPassedReads().getDuplicates());
        assertEquals(4, flagstat.qcPassedReads().getMapped());
        assertEquals(5, flagstat.qcPassedReads().getPairedInSequencing());
        assertEquals(6, flagstat.qcPassedReads().getRead1());
        assertEquals(7, flagstat.qcPassedReads().getRead2());
        assertEquals(8, flagstat.qcPassedReads().getProperlyPaired());
        assertEquals(9, flagstat.qcPassedReads().getItselfAndMateMapped());
        assertEquals(10, flagstat.qcPassedReads().getSingletons());
        assertEquals(11, flagstat.qcPassedReads().getMateMappedToDifferentChr());
        assertEquals(12, flagstat.qcPassedReads().getMateMappedToDifferentChrMapQ5());

        assertEquals(20, flagstat.qcFailedReads().getTotal());
        assertEquals(21, flagstat.qcFailedReads().getSecondary());
        assertEquals(22, flagstat.qcFailedReads().getSupplementary());
        assertEquals(23, flagstat.qcFailedReads().getDuplicates());
        assertEquals(24, flagstat.qcFailedReads().getMapped());
        assertEquals(25, flagstat.qcFailedReads().getPairedInSequencing());
        assertEquals(26, flagstat.qcFailedReads().getRead1());
        assertEquals(27, flagstat.qcFailedReads().getRead2());
        assertEquals(28, flagstat.qcFailedReads().getProperlyPaired());
        assertEquals(29, flagstat.qcFailedReads().getItselfAndMateMapped());
        assertEquals(30, flagstat.qcFailedReads().getSingletons());
        assertEquals(31, flagstat.qcFailedReads().getMateMappedToDifferentChr());
        assertEquals(32, flagstat.qcFailedReads().getMateMappedToDifferentChrMapQ5());
    }
}
