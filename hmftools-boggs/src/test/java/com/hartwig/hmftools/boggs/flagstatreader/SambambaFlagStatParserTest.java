package com.hartwig.hmftools.boggs.flagstatreader;

import com.google.common.io.Resources;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import static org.junit.Assert.assertEquals;

public class SambambaFlagStatParserTest {

    @Test
    public void canParseExampleFile() throws IOException {
        URL exampleFlagStatURL = Resources.getResource("flagstats/example.flagstat");
        String exampleFlagStatFile = exampleFlagStatURL.getPath();

        FlagStatParser parser = new SambambaFlagStatParser();
        FlagStatData flagStatData = parser.parse(new File(exampleFlagStatFile));

        assertEquals(0, flagStatData.qcPassedReads().total());
        assertEquals(1, flagStatData.qcPassedReads().secondary());
        assertEquals(2, flagStatData.qcPassedReads().supplementary());
        assertEquals(3, flagStatData.qcPassedReads().duplicates());
        assertEquals(4, flagStatData.qcPassedReads().mapped());
        assertEquals(5, flagStatData.qcPassedReads().pairedInSequencing());
        assertEquals(6, flagStatData.qcPassedReads().read1());
        assertEquals(7, flagStatData.qcPassedReads().read2());
        assertEquals(8, flagStatData.qcPassedReads().properlyPaired());
        assertEquals(9, flagStatData.qcPassedReads().itselfAndMateMapped());
        assertEquals(10, flagStatData.qcPassedReads().singletons());
        assertEquals(11, flagStatData.qcPassedReads().mateMappedToDifferentChr());
        assertEquals(12, flagStatData.qcPassedReads().mateMappedToDifferentChrMapQ5());

        assertEquals(20, flagStatData.qcFailedReads().total());
        assertEquals(21, flagStatData.qcFailedReads().secondary());
        assertEquals(22, flagStatData.qcFailedReads().supplementary());
        assertEquals(23, flagStatData.qcFailedReads().duplicates());
        assertEquals(24, flagStatData.qcFailedReads().mapped());
        assertEquals(25, flagStatData.qcFailedReads().pairedInSequencing());
        assertEquals(26, flagStatData.qcFailedReads().read1());
        assertEquals(27, flagStatData.qcFailedReads().read2());
        assertEquals(28, flagStatData.qcFailedReads().properlyPaired());
        assertEquals(29, flagStatData.qcFailedReads().itselfAndMateMapped());
        assertEquals(30, flagStatData.qcFailedReads().singletons());
        assertEquals(31, flagStatData.qcFailedReads().mateMappedToDifferentChr());
        assertEquals(32, flagStatData.qcFailedReads().mateMappedToDifferentChrMapQ5());
    }
}
