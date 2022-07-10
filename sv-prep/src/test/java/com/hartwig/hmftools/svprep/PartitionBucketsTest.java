package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.CHR_1;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.HOTSPOT_CACHE;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.READ_FILTERS;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.createSamRecord;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.readIdStr;
import static com.hartwig.hmftools.svprep.reads.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.JUNCTION;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.JunctionTracker;
import com.hartwig.hmftools.svprep.reads.PartitionBuckets;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadRecord;
import com.hartwig.hmftools.svprep.reads.BucketData;
import com.hartwig.hmftools.svprep.reads.ReadType;

import org.junit.Test;

public class PartitionBucketsTest
{
    private static final String REF_BASES = generateRandomBases(500);

    private final ChrBaseRegion mPartitionRegion;
    private final JunctionTracker mJunctionTracker;

    public PartitionBucketsTest()
    {
        mPartitionRegion = new ChrBaseRegion(CHR_1, 1, 5000);
        mJunctionTracker = new JunctionTracker(mPartitionRegion, READ_FILTERS, HOTSPOT_CACHE);
    }

    private void addRead(final ReadRecord read, final ReadType readType)
    {
        read.setReadType(readType);
        mJunctionTracker.processRead(read);
    }

    @Test
    public void testBucketJunctions()
    {
        // create junctions in the first and second buckets
        int readId = 0;

        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, REF_BASES.substring(0, 100), "30S70M"));

        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 820, REF_BASES.substring(20, 120), "100M"));

        addRead(read1, JUNCTION);
        addRead(read2, JUNCTION);

        ReadRecord suppRead1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, REF_BASES.substring(0, 73), "3S70M"));

        addRead(suppRead1, CANDIDATE_SUPPORT);

        // spanning read but with the junction still in the first bucket
        ReadRecord read3 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 950, REF_BASES.substring(0, 100), "30S70M"));

        ReadRecord read4 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 980, REF_BASES.substring(20, 120), "100M"));

        addRead(read3, JUNCTION);
        addRead(read4, JUNCTION);

        ReadRecord suppRead2 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 950, REF_BASES.substring(0, 73), "3S70M"));

        addRead(suppRead2, CANDIDATE_SUPPORT);

        // a read group starting in the first bucket but with the junction in the next
        ReadRecord read5 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 950, REF_BASES.substring(20, 120), "100M"));

        ReadRecord read6 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 980, REF_BASES.substring(0, 100), "70M30S"));

        addRead(read5, JUNCTION);
        addRead(read6, JUNCTION);

        ReadRecord suppRead3 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 990, REF_BASES.substring(0, 63), "60M3S"));

        // partitionBuckets.findBucket(suppRead3.start()).addSupportingRead(suppRead3);
        addRead(suppRead3, CANDIDATE_SUPPORT);

        // and a junction in the next bucket but with a supporting read in the previous
        ReadRecord read7 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 1010, REF_BASES.substring(10, 90), "50M30S"));

        ReadRecord read8 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 1010, REF_BASES.substring(0, 50), "50M"));

        addRead(read7, JUNCTION);
        addRead(read7, JUNCTION);

        ReadRecord suppRead4 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 990, REF_BASES.substring(0, 73), "70M3S"));

        // partitionBuckets.findBucket(readGroup1.minStartPosition()).addSupportingRead(suppRead4);
        addRead(suppRead4, CANDIDATE_SUPPORT);

        mJunctionTracker.createJunctions();

        assertEquals(4, mJunctionTracker.junctions().size());
        assertEquals(1, mJunctionTracker.junctions().get(0).supportingReadCount());
        assertEquals(1, mJunctionTracker.junctions().get(1).supportingReadCount());
        assertEquals(1, mJunctionTracker.junctions().get(2).supportingReadCount());
        assertEquals(1, mJunctionTracker.junctions().get(3).supportingReadCount());
    }

    @Test
    public void testInternalDeletes()
    {
        // initial delete is too short
        int readId = 0;

        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 10, REF_BASES.substring(0, 80), "20M10D50M"));

        addRead(read1, JUNCTION);

        // then a simple one
        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 100, REF_BASES.substring(0, 80), "20M40D20M"));

        addRead(read2, JUNCTION);

        // and a more complicated one
        // 5S10M2D10M3I10M35D10M2S from base 210: 10-19 match, 20-21 del, 22-31 match, ignore insert, 32-41 match, 42-76 del, 77-86 match

        ReadRecord read3 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 210, REF_BASES.substring(0, 1), "5S10M2D10M3I10M35D10M2S"));

        addRead(read3, JUNCTION);

        mJunctionTracker.createJunctions();

        assertEquals(4, mJunctionTracker.junctions().size());
        assertEquals(119, mJunctionTracker.junctions().get(0).Position);
        assertEquals(160, mJunctionTracker.junctions().get(1).Position);

        assertEquals(241, mJunctionTracker.junctions().get(2).Position);
        assertEquals(277, mJunctionTracker.junctions().get(3).Position);
    }
}
