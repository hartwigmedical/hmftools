package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.CHR_1;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.createSamRecord;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.readIdStr;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.junit.Test;

public class PartitionBucketsTest
{
    private static final String REF_BASES = generateRandomBases(500);

    @Test
    public void testBucketJunctions()
    {
        ChrBaseRegion partitionRegion = new ChrBaseRegion(CHR_1, 1, 5000);
        PartitionBuckets partitionBuckets = new PartitionBuckets(partitionRegion, 5000, 1000);

        // create junctions in the first and second buckets
        int readId = 1;

        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(readId++), CHR_1, 800, REF_BASES.substring(0, 100), "30S70M"));

        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 820, REF_BASES.substring(20, 120), "100M"));

        ReadGroup readGroup1 = new ReadGroup(read1);
        readGroup1.addRead(read2);

        partitionBuckets.findBucket(readGroup1.minStartPosition()).addReadGroup(readGroup1);

        // spanning read but with the junction still in the first bucket
        ReadRecord read3 = ReadRecord.from(createSamRecord(
                readIdStr(readId++), CHR_1, 950, REF_BASES.substring(0, 100), "30S70M"));

        ReadRecord read4 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 980, REF_BASES.substring(20, 120), "100M"));

        ReadGroup readGroup2 = new ReadGroup(read3);
        readGroup2.addRead(read4);

        partitionBuckets.findBucket(readGroup2.minStartPosition()).addReadGroup(readGroup2);

        // a read group starting in the first bucket but with the junction in the next
        ReadRecord read5 = ReadRecord.from(createSamRecord(
                readIdStr(readId++), CHR_1, 950, REF_BASES.substring(20, 120), "100M"));

        ReadRecord read6 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 980, REF_BASES.substring(0, 100), "70M30S"));

        ReadGroup readGroup3 = new ReadGroup(read5);
        readGroup3.addRead(read6);

        partitionBuckets.findBucket(readGroup3.minStartPosition()).addReadGroup(readGroup3);

        assertEquals(1, partitionBuckets.getBucketCount());
        SvBucket bucket1 = partitionBuckets.getBuckets().get(0);
        bucket1.createJunctions();
        partitionBuckets.transferToNext(bucket1);

        assertEquals(2, partitionBuckets.getBucketCount());

        SvBucket bucket2 = partitionBuckets.getBuckets().get(1);
        bucket2.createJunctions();
    }
}
