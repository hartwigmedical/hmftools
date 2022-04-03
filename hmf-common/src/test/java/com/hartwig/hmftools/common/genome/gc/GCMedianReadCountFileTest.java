package com.hartwig.hmftools.common.genome.gc;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GCMedianReadCountFileTest {

    private static final String GC_MEDIAN_PATH = Resources.getResource("gc/EXAMPLE.purple.gc.median").getPath();

    @Test
    public void canLoadWithoutExtendingRegion() throws IOException {
        final GCMedianReadCount readCount = GCMedianReadCountFile.read(false, GC_MEDIAN_PATH);
        testMinMax(readCount, -1, -1);
    }

    @Test
    public void canLoadAndExtendRegion() throws IOException {
        final GCMedianReadCount readCount = GCMedianReadCountFile.read(true, GC_MEDIAN_PATH);
        testMinMax(readCount, 1144, 880);
    }

    private static void testMinMax(@NotNull final GCMedianReadCount victim, int expectedMin, int expectedMax) {
        for (int i = 0; i < 20; i++) {
            assertEquals(expectedMin, victim.medianReadCount(new ImmutableGCBucket(i)), 1e-10);
        }

        assertEquals(1144, victim.medianReadCount(new ImmutableGCBucket(20)), 1e-10);
        assertEquals(994, victim.medianReadCount(new ImmutableGCBucket(30)), 1e-10);
        assertEquals(-1, victim.medianReadCount(new ImmutableGCBucket(40)), 1e-10);
        assertEquals(894, victim.medianReadCount(new ImmutableGCBucket(50)), 1e-10);
        assertEquals(880, victim.medianReadCount(new ImmutableGCBucket(60)), 1e-10);

        for (int i = 61; i <= 100; i++) {
            assertEquals(expectedMax, victim.medianReadCount(new ImmutableGCBucket(i)), 1e-10);
        }
    }
}
