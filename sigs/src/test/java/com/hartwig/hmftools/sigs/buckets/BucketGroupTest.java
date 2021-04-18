package com.hartwig.hmftools.sigs.buckets;

import static junit.framework.TestCase.assertEquals;

import org.junit.Test;

public class BucketGroupTest
{

    @Test
    public void testSampleAddAndRemove()
    {
        BucketGroup bucketGroup = new BucketGroup(0);

        bucketGroup.addBucket(0, true);
        bucketGroup.addBucket(1, true);
        assertEquals(2, bucketGroup.getBucketCount());

        SampleData sample1 = new SampleData(0);
        double[] counts = {10, 10};
        sample1.setBucketCounts(counts);

        SampleData sample2 = new SampleData(1);
        double[] counts2 = {20, 20};
        sample2.setBucketCounts(counts2);

        bucketGroup.addSample(sample1.Id, counts);
        bucketGroup.addSample(sample2.Id, counts2);

        assertEquals(2, bucketGroup.getSampleCount());
        assertEquals(60.0, bucketGroup.getTotalCount());
        assertEquals(0, bucketGroup.getSampleIndex(sample1.Id));
        assertEquals(1, bucketGroup.getSampleIndex(sample2.Id));

        assertEquals(20.0, bucketGroup.getSampleCount(sample1.Id));
        assertEquals(40.0, bucketGroup.getSampleCount(sample2.Id));
        assertEquals(30.0, bucketGroup.getBucketCounts()[0]);
        assertEquals(30.0, bucketGroup.getBucketCounts()[1]);

        // add more for the first sample
        bucketGroup.addSampleCounts(0, counts);
        assertEquals(80.0, bucketGroup.getTotalCount());
        assertEquals(40.0, bucketGroup.getSampleCount(sample1.Id));
        assertEquals(40.0, bucketGroup.getBucketCounts()[0]);
        assertEquals(40.0, bucketGroup.getBucketCounts()[1]);

        bucketGroup.removeSampleAllocation(sample1, 0, false);
        assertEquals(1, bucketGroup.getSampleCount());
        assertEquals(40.0, bucketGroup.getTotalCount());
        assertEquals(0, bucketGroup.getSampleIndex(sample2.Id));
        assertEquals(20.0, bucketGroup.getBucketCounts()[0]);
        assertEquals(20.0, bucketGroup.getBucketCounts()[1]);
    }

}
