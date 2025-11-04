package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.calculations.GCPail.bucketIndex;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.cobalt.count.DepthReading;

import org.junit.Test;
import org.mockito.Mockito;

public class CobaltWindowTest
{
    DepthReading readDepth = new DepthReading("1", 1001, 82, 0.49);
    DepthReading readDepth2 = new DepthReading("2", 5001, 1892, 0.52);
    DepthReading readDepthY = new DepthReading("Y", 5001, 1892, 0.52);

    @Test
    public void construct()
    {
        CobaltWindow window = new CobaltWindow(_1, readDepth, true, true);
        assertEquals(readDepth.StartPosition, window.Position);
        assertTrue(window.IsInExcludedRegion);
        assertTrue(window.IsInTargetRegion);
        assertNull(window.GcBucket);
    }

    @Test
    public void constructWithBucket()
    {
        GCPail bucket = new GCPail(49);
        CobaltWindow window = new CobaltWindow(_1, readDepth.StartPosition, readDepth.ReadDepth, readDepth.ReadGcContent, bucket, true);
        assertEquals(readDepth.StartPosition, window.Position);
        assertFalse(window.IsInExcludedRegion);
        assertTrue(window.IsInTargetRegion);
        assertEquals(bucket, window.GcBucket);
    }

    @Test
    public void included()
    {
        GCPail bucket = new GCPail(49);
        CobaltWindow window = new CobaltWindow(_1, readDepth.StartPosition, readDepth.ReadDepth, readDepth.ReadGcContent, bucket, true);
        assertTrue(window.include());

        window = new CobaltWindow(_2, readDepth2, true, true);
        assertFalse(window.include());
        window = new CobaltWindow(_2, readDepth2, true, false);
        assertFalse(window.include());
        window = new CobaltWindow(_2, readDepth2, false, true);
        assertTrue(window.include());
        window = new CobaltWindow(_2, readDepth2, false, false);
        assertFalse(window.include());
    }

    @Test
    public void toBamRatio()
    {
        GCPail bucket = new GCPail(49);
        CobaltWindow window = new CobaltWindow(_1, readDepth.StartPosition, readDepth.ReadDepth, readDepth.ReadGcContent, bucket, true);
        BamRatio ratio = window.toBamRatio();
        assertEquals(readDepth.StartPosition, window.Position);
        assertEquals(readDepth.ReadDepth, ratio.readDepth(), 0.001);
        assertEquals(readDepth.ReadGcContent, ratio.gcContent(), 0.001);

        window = new CobaltWindow(_2, readDepth2, true, true);
        ratio = window.toBamRatio();
        assertEquals(readDepth2.StartPosition, window.Position);
        assertEquals(readDepth2.ReadDepth, ratio.readDepth(), 0.001);
        assertEquals(-1.0, ratio.ratio(), 0.001); // excluded
        assertEquals(readDepth2.ReadGcContent, ratio.gcContent(), 0.001);

        window = new CobaltWindow(_2, readDepth2, true, false);
        ratio = window.toBamRatio();
        assertEquals(readDepth2.StartPosition, window.Position);
        assertEquals(readDepth2.ReadDepth, ratio.readDepth(), 0.001);
        assertEquals(-1.0, ratio.ratio(), 0.001);
        assertEquals(readDepth2.ReadGcContent, ratio.gcContent(), 0.001);

        window = new CobaltWindow(_2, readDepth2, false, true);
        ratio = window.toBamRatio();
        assertEquals(readDepth2.StartPosition, window.Position);
        assertEquals(readDepth2.ReadDepth, ratio.readDepth(), 0.001);
        assertEquals(readDepth2.ReadDepth, ratio.ratio(), 0.001);
        assertEquals(readDepth2.ReadGcContent, ratio.gcContent(), 0.001);

        window = new CobaltWindow(_2, readDepth2, false, false);
        ratio = window.toBamRatio();
        assertEquals(readDepth2.StartPosition, window.Position);
        assertEquals(readDepth2.ReadDepth, ratio.readDepth(), 0.001);
        assertEquals(-1.0, ratio.ratio(), 0.001);
        assertEquals(readDepth2.ReadGcContent, ratio.gcContent(), 0.001);
    }

    @Test
    public void correctedByReferenceValueTest()
    {
        WindowStatuses filter = Mockito.mock(WindowStatuses.class);
        Mockito.when(filter.referenceGcValueForWindow(_1, readDepth.StartPosition)).thenReturn(0.89);

        CobaltWindow window = new CobaltWindow(_1, readDepth.StartPosition, readDepth.ReadDepth, readDepth.ReadGcContent, null, true);
        CobaltWindow corrected = window.correctedByReferenceValue(filter);
        Mockito.verify(filter, Mockito.never()).referenceGcValueForWindow(_1, readDepth.StartPosition);
        assertSame(corrected, window);

        window = new CobaltWindow(_1, readDepth.StartPosition, -1.0, readDepth.ReadGcContent, null, true);
        corrected = window.correctedByReferenceValue(filter);
        Mockito.verify(filter, Mockito.never()).referenceGcValueForWindow(_1, readDepth.StartPosition);
        assertSame(corrected, window);

        CobaltWindow window2 = new CobaltWindow(_1, readDepth.StartPosition, 0.99, readDepth.ReadGcContent, null, true);
        CobaltWindow corrected2 = window2.correctedByReferenceValue(filter);
        Mockito.verify(filter, Mockito.only()).referenceGcValueForWindow(_1, readDepth.StartPosition);
        assertEquals(0.89, corrected2.GcContent, 0.001);
    }

    @Test
    public void bucketedTest()
    {
        final GCPail gcBucket = new GCPail(bucketIndex(readDepth2.ReadGcContent));
        assertEquals(0.0, gcBucket.median(), 0.001);

        CobaltWindow window = new CobaltWindow(_2, readDepth2, false, true);
        CobaltWindow bucketed = window.bucketed(gcBucket);
        assertEquals(gcBucket, bucketed.GcBucket);
        assertFalse(bucketed.IsInExcludedRegion);
        assertTrue(bucketed.IsInTargetRegion);
        assertEquals(readDepth2.ReadDepth, gcBucket.median(), 0.001);
    }

    @Test
    public void bucketedAllosomeTest()
    {
        final GCPail gcBucket = new GCPail(bucketIndex(readDepthY.ReadGcContent));
        assertEquals(0.0, gcBucket.median(), 0.001);

        CobaltWindow window = new CobaltWindow(_Y, readDepthY, false, true);
        CobaltWindow bucketed = window.bucketed(gcBucket);
        assertEquals(gcBucket, bucketed.GcBucket);
        assertFalse(bucketed.IsInExcludedRegion);
        assertTrue(bucketed.IsInTargetRegion);
        assertEquals(0.0, gcBucket.median(), 0.001);
    }

    @Test
    public void bucketedOffTargetRegionTest()
    {
        final GCPail gcBucket = new GCPail(bucketIndex(readDepth2.ReadGcContent));
        assertEquals(0.0, gcBucket.median(), 0.001);

        CobaltWindow window = new CobaltWindow(_2, readDepth2, false, false);
        CobaltWindow bucketed = window.bucketed(gcBucket);
        assertEquals(gcBucket, bucketed.GcBucket);
        assertFalse(bucketed.IsInExcludedRegion);
        assertFalse(bucketed.IsInTargetRegion);
        assertEquals(0.0, gcBucket.median(), 0.001);
    }

    @Test
    public void bucketedDepth0Test()
    {
        final GCPail gcBucket = new GCPail(bucketIndex(0.0));
        gcBucket.addReading(99.0);
        assertEquals(99.0, gcBucket.median(), 0.001);

        DepthReading depth0 = new DepthReading("2", 5001, 0.0, Double.NaN);
        CobaltWindow window = new CobaltWindow(_2, depth0, false, true);
        CobaltWindow bucketed = window.bucketed(gcBucket);
        assertEquals(gcBucket, bucketed.GcBucket);
        assertFalse(bucketed.IsInExcludedRegion);
        assertTrue(bucketed.IsInTargetRegion);
        assertEquals(99.0, gcBucket.median(), 0.001);
    }

    @Test
    public void testEquals()
    {
        CobaltWindow window1 = new CobaltWindow(_1, readDepth, true, true);
        CobaltWindow window2 = new CobaltWindow(_1, readDepth, false, false);
        CobaltWindow window3 = new CobaltWindow(_1, readDepth, false, true);
        assertEquals(window1, window2);
        assertEquals(window1, window3);
        assertEquals(window2, window3);

        DepthReading differentPosition = new DepthReading("1", 2001, 82, 0.49);
        CobaltWindow window4 = new CobaltWindow(_1, differentPosition, true, true);
        assertNotEquals(window1, window4);
        CobaltWindow window5 = new CobaltWindow(_2, readDepth2, true, true);
        assertNotEquals(window1, window5);
        assertNotEquals(null, window1);
        assertNotEquals("Not a CobaltWindow", window1);
    }

    @Test
    public void testHashCode()
    {
        CobaltWindow window1 = new CobaltWindow(_1, readDepth, true, true);
        CobaltWindow window2 = new CobaltWindow(_1, readDepth, false, false);

        assertEquals(window1.hashCode(), window2.hashCode());
    }
}
