package com.hartwig.hmftools.svanalysis.visualisation.circos;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CircosDataWriterTest {

    @Test
    public void testShorthand() {
        assertEquals("1", CircosDataWriter.shorthand(1));
        assertEquals("99", CircosDataWriter.shorthand(99));
        assertEquals("0.1k", CircosDataWriter.shorthand(100));
        assertEquals("1.0k", CircosDataWriter.shorthand(1_000));
        assertEquals("99.9k", CircosDataWriter.shorthand(99_949));
        assertEquals("0.1m", CircosDataWriter.shorthand(99_950));
        assertEquals("0.1m", CircosDataWriter.shorthand(100_000));
    }

}
