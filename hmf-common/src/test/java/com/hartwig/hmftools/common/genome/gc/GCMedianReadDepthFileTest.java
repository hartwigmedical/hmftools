package com.hartwig.hmftools.common.genome.gc;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class GCMedianReadDepthFileTest
{
    private static final double DELTA = 1e-2;

    @Test
    public void testCanLoad() throws IOException
    {
        @SuppressWarnings("UnstableApiUsage")
        final String GC_MEDIAN_PATH = Resources.getResource("gc/EXAMPLE.cobalt.gc.median").getPath();

        final GCMedianReadDepth readDepth = GCMedianReadDepthFile.read(GC_MEDIAN_PATH);
        validateReadDepth(readDepth);
    }

    @Test
    public void testCanLoadOldVersion() throws IOException
    {
        @SuppressWarnings("UnstableApiUsage")
        final String GC_MEDIAN_PATH = Resources.getResource("gc/EXAMPLE.cobalt.gc.median.oldversion").getPath();

        final GCMedianReadDepth readDepth = GCMedianReadDepthFile.read(GC_MEDIAN_PATH);
        validateReadDepth(readDepth);
    }

    private void validateReadDepth(final GCMedianReadDepth readDepth)
    {
        assertEquals(145.11, readDepth.meanReadDepth(), DELTA);
        assertEquals(141.79, readDepth.medianReadDepth(), DELTA);

        assertEquals(150.09, readDepth.medianReadDepth(new ImmutableGCBucket(30)), DELTA);
        assertEquals(165.5, readDepth.medianReadDepth(new ImmutableGCBucket(22)), DELTA);
        assertEquals(150.09, readDepth.medianReadDepth(new ImmutableGCBucket(30)), DELTA);
        assertEquals(-1, readDepth.medianReadDepth(new ImmutableGCBucket(40)), DELTA);
        assertEquals(134.99, readDepth.medianReadDepth(new ImmutableGCBucket(50)), DELTA);
        assertEquals(132.88, readDepth.medianReadDepth(new ImmutableGCBucket(60)), DELTA);
    }
}
