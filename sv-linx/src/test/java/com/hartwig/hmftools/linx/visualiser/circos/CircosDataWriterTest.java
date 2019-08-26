package com.hartwig.hmftools.linx.visualiser.circos;

import static com.hartwig.hmftools.linx.visualiser.circos.CircosDataWriter.shorthand;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CircosDataWriterTest
{

    @Test
    public void testShorthand()
    {
        assertEquals("1", shorthand(1));
        assertEquals("99", shorthand(99));
        assertEquals("0.1k", shorthand(100));
        assertEquals("1.0k", shorthand(1_000));
        assertEquals("99.9k", shorthand(99_949));
        assertEquals("0.1m", shorthand(99_950));
        assertEquals("0.1m", shorthand(100_000));
    }

//    @Test
//    public void testThickness() {
//        assertEquals(1, thicknessPixels(0.5, 6), 0.01);
//        assertEquals(2, thicknessPixels(1, 6), 0.01);
//        assertEquals(4, thicknessPixels(2, 6), 0.01);
//        assertEquals(6, thicknessPixels(3, 6), 0.01);
//        assertEquals(8, thicknessPixels(4, 6), 0.01);
//        assertEquals(10, thicknessPixels(5, 6), 0.01);
//        assertEquals(12, thicknessPixels(6, 6), 0.01);
//
//
//        assertEquals(1, thicknessPixels(0.5, 12), 0.01);
//        assertEquals(1, thicknessPixels(1, 12), 0.01);
//        assertEquals(2, thicknessPixels(2, 12), 0.01);
//        assertEquals(3, thicknessPixels(3, 12), 0.01);
//        assertEquals(4, thicknessPixels(4, 12), 0.01);
//        assertEquals(5, thicknessPixels(5, 12), 0.01);
//        assertEquals(6, thicknessPixels(6, 12), 0.01);
//        assertEquals(7, thicknessPixels(7, 12), 0.01);
//        assertEquals(12, thicknessPixels(12, 12), 0.01);
//        assertEquals(12, thicknessPixels(100, 12), 0.01);
//    }
}
