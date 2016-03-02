package com.hartwig.hmftools.sullivan;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class DefaultNameConverterTest {

    @Test
    public void converterWorksAsExpected() {
        String originalFileName = "CPCT02010015R_BHKCWVCCXX_S3_L001_R1_001.fastq.gz";
        String expectedConvertedName = "CPCT02010015R_BHKCWVCCXX_S3_L001_001_1.fastq";

        FileNameConverter converter = new DefaultNameConverter();

        assertEquals(expectedConvertedName, converter.apply(originalFileName));
    }
}
