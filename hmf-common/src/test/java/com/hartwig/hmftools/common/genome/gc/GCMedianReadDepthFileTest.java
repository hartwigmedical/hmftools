package com.hartwig.hmftools.common.genome.gc;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GCMedianReadDepthFileTest
{
    @SuppressWarnings("UnstableApiUsage")
    private static final String GC_MEDIAN_PATH = Resources.getResource("gc/EXAMPLE.purple.gc.median").getPath();

    @Test
    public void canLoadWithoutExtendingRegion() throws IOException
    {
        final GCMedianReadDepth readCount = GCMedianReadDepthFile.read(false, GC_MEDIAN_PATH);
        testMinMax(readCount, -1, -1);
    }

    @Test
    public void canLoadAndExtendRegion() throws IOException
    {
        final GCMedianReadDepth readCount = GCMedianReadDepthFile.read(true, GC_MEDIAN_PATH);
        testMinMax(readCount, 1144, 880);
    }

    private static void testMinMax(@NotNull final GCMedianReadDepth victim, int expectedMin, int expectedMax)
    {
        for(int i = 0; i < 20; i++)
        {
            assertEquals(expectedMin, victim.medianReadDepth(new ImmutableGCBucket(i)), 1e-10);
        }

        assertEquals(1144, victim.medianReadDepth(new ImmutableGCBucket(20)), 1e-10);
        assertEquals(994, victim.medianReadDepth(new ImmutableGCBucket(30)), 1e-10);
        assertEquals(-1, victim.medianReadDepth(new ImmutableGCBucket(40)), 1e-10);
        assertEquals(894, victim.medianReadDepth(new ImmutableGCBucket(50)), 1e-10);
        assertEquals(880, victim.medianReadDepth(new ImmutableGCBucket(60)), 1e-10);

        for(int i = 61; i <= 100; i++)
        {
            assertEquals(expectedMax, victim.medianReadDepth(new ImmutableGCBucket(i)), 1e-10);
        }
    }
}
