package com.hartwig.hmftools.boggs.flagstatreader;

import com.google.common.io.Resources;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import static org.junit.Assert.assertEquals;

public class SambambaFlagstatParserTest2 {

    @Test
    public void canParseExampleFile() throws IOException {
        URL exampleFlagstatURL = Resources.getResource("flagstats/example.flagstat");
        String exampleFlagstatFile = exampleFlagstatURL.getPath();

        FlagStatParser parser = new SambambaFlagStatParser();
        FlagStatData flagstatData = parser.parse(new File(exampleFlagstatFile));

        assertEquals(0, flagstatData.qcPassedReads().total());
        assertEquals(1, flagstatData.qcPassedReads().secondary());
        assertEquals(2, flagstatData.qcPassedReads().supplementary());
        assertEquals(3, flagstatData.qcPassedReads().duplicates());
        assertEquals(4, flagstatData.qcPassedReads().mapped());
        assertEquals(5, flagstatData.qcPassedReads().pairedInSequencing());
        assertEquals(6, flagstatData.qcPassedReads().read1());
        assertEquals(7, flagstatData.qcPassedReads().read2());
        assertEquals(8, flagstatData.qcPassedReads().properlyPaired());
        assertEquals(9, flagstatData.qcPassedReads().itselfAndMateMapped());
        assertEquals(10, flagstatData.qcPassedReads().singletons());
        assertEquals(11, flagstatData.qcPassedReads().mateMappedToDifferentChr());
        assertEquals(12, flagstatData.qcPassedReads().mateMappedToDifferentChrMapQ5());

        assertEquals(20, flagstatData.qcFailedReads().total());
        assertEquals(21, flagstatData.qcFailedReads().secondary());
        assertEquals(22, flagstatData.qcFailedReads().supplementary());
        assertEquals(23, flagstatData.qcFailedReads().duplicates());
        assertEquals(24, flagstatData.qcFailedReads().mapped());
        assertEquals(25, flagstatData.qcFailedReads().pairedInSequencing());
        assertEquals(26, flagstatData.qcFailedReads().read1());
        assertEquals(27, flagstatData.qcFailedReads().read2());
        assertEquals(28, flagstatData.qcFailedReads().properlyPaired());
        assertEquals(29, flagstatData.qcFailedReads().itselfAndMateMapped());
        assertEquals(30, flagstatData.qcFailedReads().singletons());
        assertEquals(31, flagstatData.qcFailedReads().mateMappedToDifferentChr());
        assertEquals(32, flagstatData.qcFailedReads().mateMappedToDifferentChrMapQ5());
    }
}
