package com.hartwig.hmftools.sullivan;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class AncientNameConverterTest {

    @Test
    public void convertsAsExpected() {
        String originalFileName = "MMC01013-00A00A0_S4_L001_R1_001.fastq.gz";
        String expectedConvertedName = "MMC01013-00A00A0_S4_L001_001_1.fastq";

        FileNameConverter converter = new AncientNameConverter();

        assertEquals(expectedConvertedName, converter.apply(originalFileName));
    }
}
