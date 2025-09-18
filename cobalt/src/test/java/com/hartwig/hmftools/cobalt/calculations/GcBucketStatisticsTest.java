package com.hartwig.hmftools.cobalt.calculations;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class GcBucketStatisticsTest
{
    @Test
    public void mediansWithNoDataTest()
    {
        GCPailsList pailsList = new GCPailsList();
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        for (int i=0; i<100; i++)
        {
            assertEquals(-1.0, bucketStatistics.medianReadDepth(i), 0.001);
        }
    }

    @Test
    public void mediansOneBucketOneReadingTest()
    {
        GCPailsList pailsList = new GCPailsList();
        pailsList.getGCPail(0.5).addReading(10);
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        assertEquals(10, bucketStatistics.medianReadDepth(50), 0.001);
    }

    @Test
    public void mediansOneBucketMultipleReadingsTest()
    {
        GCPailsList pailsList = new GCPailsList();
        pailsList.getGCPail(0.5).addReading(10);
        pailsList.getGCPail(0.5).addReading(10);
        pailsList.getGCPail(0.5).addReading(12);
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        assertEquals(10, bucketStatistics.medianReadDepth(50), 0.001);
    }

    @Test
    public void mediansTwoAdjoiningBucketTest()
    {
        GCPailsList pailsList = new GCPailsList();
        pailsList.getGCPail(0.50).addReading(8.0);
        pailsList.getGCPail(0.51).addReading(9.0);
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        assertEquals(0.0, bucketStatistics.medianReadDepth(49), 0.001);
        assertEquals(8.5, bucketStatistics.medianReadDepth(50), 0.001);
        assertEquals(8.5, bucketStatistics.medianReadDepth(51), 0.001);
        assertEquals(0.0, bucketStatistics.medianReadDepth(52), 0.001);
    }

    @Test
    public void mediansTwoBucketsWithAGapTest()
    {
        GCPailsList pailsList = new GCPailsList();
        pailsList.getGCPail(0.50).addReading(8.0);
        pailsList.getGCPail(0.52).addReading(10.0);
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        assertEquals(0.0, bucketStatistics.medianReadDepth(49), 0.001);
        assertEquals(8.0, bucketStatistics.medianReadDepth(50), 0.001);
        assertEquals(0.0, bucketStatistics.medianReadDepth(51), 0.001);
        assertEquals(10.0, bucketStatistics.medianReadDepth(52), 0.001);
        assertEquals(0.0, bucketStatistics.medianReadDepth(53), 0.001);
    }

    @Test
    public void manyBucketsTest()
    {
        GCPailsList pailsList = new GCPailsList();
        // Median for 50 is 8
        pailsList.getGCPail(0.50).addReading(8.0);
        pailsList.getGCPail(0.50).addReading(8.0);
        pailsList.getGCPail(0.50).addReading(10.0);

        // median for 51 is 10
        pailsList.getGCPail(0.51).addReading(12.0);
        pailsList.getGCPail(0.51).addReading(8.0);

        // median for 52 is 12
        pailsList.getGCPail(0.52).addReading(12.0);
        pailsList.getGCPail(0.52).addReading(12.0);

        // median for 53 is 14
        pailsList.getGCPail(0.53).addReading(15.0);
        pailsList.getGCPail(0.53).addReading(14.0);
        pailsList.getGCPail(0.53).addReading(13.0);

        // median for 54 is 16
        pailsList.getGCPail(0.54).addReading(16.0);

        // median for 55 is 18
        pailsList.getGCPail(0.55).addReading(18.0);

        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        assertEquals(0.0, bucketStatistics.medianReadDepth(49), 0.001);
        assertEquals(9.0, bucketStatistics.medianReadDepth(50), 0.001);
        assertEquals(10.0, bucketStatistics.medianReadDepth(51), 0.001);
        assertEquals(12.0, bucketStatistics.medianReadDepth(52), 0.001);
        assertEquals(14.0, bucketStatistics.medianReadDepth(53), 0.001);
        assertEquals(16.0, bucketStatistics.medianReadDepth(54), 0.001);
        assertEquals(17.0, bucketStatistics.medianReadDepth(55), 0.001);
        assertEquals(0.0, bucketStatistics.medianReadDepth(56), 0.001);
    }

    @Test
    public void allBucketsTest()
    {
        GCPailsList pailsList = new GCPailsList();
        for (int i=0; i<100; i++)
        {
            pailsList.getGCPail(0.01 * i).addReading(1.0 * i);
        }

        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 20, 80);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(0), 0.001);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(10), 0.001);
        assertEquals(19.0, bucketStatistics.medianReadDepth(19), 0.001);
        assertEquals(20.0, bucketStatistics.medianReadDepth(20), 0.001);
        assertEquals(50.0, bucketStatistics.medianReadDepth(50), 0.001);
        assertEquals(79.0, bucketStatistics.medianReadDepth(79), 0.001);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(80), 0.001);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(81), 0.001);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(99), 0.001);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(100), 0.001);
    }

    @Test
    public void isAllowedTest()
    {
        GCPailsList pailsList = new GCPailsList();
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        for (int i=0; i<100; i++)
        {
            if (i < 23 || i > 67)
            {
                assertFalse(bucketStatistics.isAllowed(new GCPail(i)));
            }
            else
            {
                assertTrue(bucketStatistics.isAllowed(new GCPail(i)));
            }
        }
    }
}
