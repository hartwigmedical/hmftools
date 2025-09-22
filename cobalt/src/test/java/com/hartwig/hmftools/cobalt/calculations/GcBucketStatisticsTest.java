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
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 20, 70);
        for (int i=0; i<100; i++)
        {
            double expectedMedian = (i < 20 || i > 69) ? -1.0 : 0.0;
            assertEquals(expectedMedian, bucketStatistics.medianReadDepth(new GCPail(i)), 0.001);
        }
    }

    @Test
    public void mediansOneBucketOneReadingTest()
    {
        GCPailsList pailsList = new GCPailsList();
        pailsList.getGCPail(0.5).addReading(10);
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        assertEquals(0.0, bucketStatistics.medianReadDepth(new GCPail(48)), 0.001);
        assertEquals(3.3333, bucketStatistics.medianReadDepth(new GCPail(49)), 0.001);
        assertEquals(3.3333, bucketStatistics.medianReadDepth(new GCPail(50)), 0.001);
        assertEquals(3.3333, bucketStatistics.medianReadDepth(new GCPail(51)), 0.001);
        assertEquals(0.0, bucketStatistics.medianReadDepth(new GCPail(52)), 0.001);
    }

    @Test
    public void mediansOneBucketMultipleReadingsTest()
    {
        GCPailsList pailsList = new GCPailsList();
        pailsList.getGCPail(0.5).addReading(10);
        pailsList.getGCPail(0.5).addReading(10);
        pailsList.getGCPail(0.5).addReading(12);
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        assertEquals(0.0, bucketStatistics.medianReadDepth(new GCPail(48)), 0.001);
        assertEquals(3.3333, bucketStatistics.medianReadDepth(new GCPail(49)), 0.001);
        assertEquals(3.3333, bucketStatistics.medianReadDepth(new GCPail(50)), 0.001);
        assertEquals(3.3333, bucketStatistics.medianReadDepth(new GCPail(51)), 0.001);
        assertEquals(0.0, bucketStatistics.medianReadDepth(new GCPail(52)), 0.001);
    }

    @Test
    public void mediansTwoAdjoiningBucketTest()
    {
        GCPailsList pailsList = new GCPailsList();
        pailsList.getGCPail(0.50).addReading(8.0);
        pailsList.getGCPail(0.51).addReading(9.0);
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        assertEquals(0.0, bucketStatistics.medianReadDepth(new GCPail(48)), 0.001);
        assertEquals(8.0/3.0, bucketStatistics.medianReadDepth(new GCPail(49)), 0.001);
        assertEquals(17.0/3.0, bucketStatistics.medianReadDepth(new GCPail(50)), 0.001);
        assertEquals(17.0/3.0, bucketStatistics.medianReadDepth(new GCPail(51)), 0.001);
        assertEquals(3.0, bucketStatistics.medianReadDepth(new GCPail(52)), 0.001);
        assertEquals(0.0, bucketStatistics.medianReadDepth(new GCPail(53)), 0.001);
    }

    @Test
    public void mediansTwoBucketsWithAGapTest()
    {
        GCPailsList pailsList = new GCPailsList();
        pailsList.getGCPail(0.50).addReading(8.0);
        pailsList.getGCPail(0.52).addReading(10.0);
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 23, 67);
        assertEquals(0.0, bucketStatistics.medianReadDepth(new GCPail(48)), 0.001);
        assertEquals(8.0/3.0, bucketStatistics.medianReadDepth(new GCPail(49)), 0.001);
        assertEquals(8.0/3.0, bucketStatistics.medianReadDepth(new GCPail(50)), 0.001);
        assertEquals(6.0, bucketStatistics.medianReadDepth(new GCPail(51)), 0.001);
        assertEquals(10.0/3.0, bucketStatistics.medianReadDepth(new GCPail(52)), 0.001);
        assertEquals(10.0/3.0, bucketStatistics.medianReadDepth(new GCPail(53)), 0.001);
        assertEquals(0.0, bucketStatistics.medianReadDepth(new GCPail(54)), 0.001);
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
        assertEquals(0.0, bucketStatistics.medianReadDepth(new GCPail(48)), 0.001);
        assertEquals(8.0/3.0, bucketStatistics.medianReadDepth(new GCPail(49)), 0.001);
        assertEquals(6.0, bucketStatistics.medianReadDepth(new GCPail(50)), 0.001);
        assertEquals(10.0, bucketStatistics.medianReadDepth(new GCPail(51)), 0.001);
        assertEquals(12.0, bucketStatistics.medianReadDepth(new GCPail(52)), 0.001);
        assertEquals(14.0, bucketStatistics.medianReadDepth(new GCPail(53)), 0.001);
        assertEquals(16.0, bucketStatistics.medianReadDepth(new GCPail(54)), 0.001);
        assertEquals(34.0/3.0, bucketStatistics.medianReadDepth(new GCPail(55)), 0.001);
        assertEquals(6.0, bucketStatistics.medianReadDepth(new GCPail(56)), 0.001);
        assertEquals(0.0, bucketStatistics.medianReadDepth(new GCPail(57)), 0.001);
    }

    @Test
    public void allBucketsTest()
    {
        GCPailsList pailsList = new GCPailsList();
        for (int i=0; i<100; i++)
        {
            pailsList.getGCPail(0.01 * i).addReading(1.0 * i);
        }

        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 19, 80);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(new GCPail(0)), 0.001);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(new GCPail(10)), 0.001);
        assertEquals(19.0, bucketStatistics.medianReadDepth(new GCPail(19)), 0.001);
        assertEquals(20.0, bucketStatistics.medianReadDepth(new GCPail(20)), 0.001);
        assertEquals(50.0, bucketStatistics.medianReadDepth(new GCPail(50)), 0.001);
        assertEquals(79.0, bucketStatistics.medianReadDepth(new GCPail(79)), 0.001);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(new GCPail(80)), 0.001);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(new GCPail(81)), 0.001);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(new GCPail(99)), 0.001);
        assertEquals(-1.0, bucketStatistics.medianReadDepth(new GCPail(100)), 0.001);
    }

    @Test
    public void isAllowedTest()
    {
        GCPailsList pailsList = new GCPailsList();
        GcBucketStatistics bucketStatistics = new GcBucketStatistics(pailsList, 24, 68);
        for (int i=0; i<100; i++)
        {
            if (i < 24 || i > 68)
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
