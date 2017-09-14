package com.hartwig.hmftools.common.gc;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GCMedianReadCountFileTest {

    private static final String BASE_PATH = Resources.getResource("gc").getPath() + File.separator;

    @Test
    public void canLoadWithouExtendingRegion() throws IOException, HartwigException {
        final GCMedianReadCount readCount = GCMedianReadCountFile.read(false, BASE_PATH + "EXAMPLE.purple.gc.median");
        testMinMax(readCount, -1, -1);
    }

    @Test
    public void canLoadAndExtendRegion() throws IOException, HartwigException {
        final GCMedianReadCount readCount = GCMedianReadCountFile.read(true, BASE_PATH + "EXAMPLE.purple.gc.median");
        testMinMax(readCount, 1144, 880);
    }

    private static void testMinMax(@NotNull final GCMedianReadCount victim, int expectedMin, int expectedMax) throws IOException, HartwigException {

        for (int i = 0; i < 20; i++) {
            assertEquals(expectedMin, victim.medianReadCount(new ImmutableGCBucket(i)));
        }

        assertEquals(1144, victim.medianReadCount(new ImmutableGCBucket(20)));
        assertEquals(994, victim.medianReadCount(new ImmutableGCBucket(30)));
        assertEquals(-1, victim.medianReadCount(new ImmutableGCBucket(40)));
        assertEquals(894, victim.medianReadCount(new ImmutableGCBucket(50)));
        assertEquals(880, victim.medianReadCount(new ImmutableGCBucket(60)));

        for (int i = 61; i <= 100; i++) {
            assertEquals(expectedMax, victim.medianReadCount(new ImmutableGCBucket(i)));
        }
    }

}
