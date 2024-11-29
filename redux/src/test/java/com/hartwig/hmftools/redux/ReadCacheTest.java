package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.redux.common.FragmentCoordReads;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadCacheTest
{
    private static final ReadCache mReadCache = new ReadCache(100, 100,false);

    @Test
    public void testBasics()
    {
        SAMRecord read1 = createRecord(CHR_1, 50, false);

        mReadCache.processRead(read1);

        assertEquals(1, mReadCache.cachedReadGroups());
        assertEquals(1, mReadCache.cachedFragCoordGroups());
        assertEquals(1, mReadCache.cachedReadCount());

        FragmentCoordReads fragCoordReads = mReadCache.popReads();
        assertNull(fragCoordReads);

        SAMRecord read2a = createRecord(CHR_1, 80, false);
        SAMRecord read2b = createRecord(CHR_1, 80, false);

        mReadCache.processRead(read2a);
        mReadCache.processRead(read2b);

        SAMRecord read3 = createRecord(CHR_1, 120, false);

        mReadCache.processRead(read3);

        assertEquals(2, mReadCache.cachedReadGroups());
        assertEquals(4, mReadCache.cachedReadCount());
        assertEquals(3, mReadCache.cachedFragCoordGroups());

        fragCoordReads = mReadCache.popReads();
        assertNull(fragCoordReads);

        SAMRecord read4 = createRecord(CHR_1, 199, false);

        mReadCache.processRead(read4);

        assertEquals(2, mReadCache.cachedReadGroups());

        // the next read triggers the first group to be processed
        SAMRecord read5 = createRecord(CHR_1, 201, false);

        mReadCache.processRead(read5);

        fragCoordReads = mReadCache.popReads();
        assertNotNull(fragCoordReads);

        assertEquals(1, fragCoordReads.SingleReads.size());
        assertEquals(1, fragCoordReads.DuplicateGroups.size());

        assertEquals(2, mReadCache.cachedReadGroups());
        assertEquals(3, mReadCache.cachedReadCount());
        assertEquals(3, mReadCache.cachedFragCoordGroups());

        // next read is on a new chromosome so clears all previous
        SAMRecord read6 = createRecord(CHR_2, 50, false);

        mReadCache.processRead(read6);

        fragCoordReads = mReadCache.popReads();
        assertNotNull(fragCoordReads);

        assertEquals(3, fragCoordReads.SingleReads.size());
        assertEquals(0, fragCoordReads.DuplicateGroups.size());

        assertEquals(1, mReadCache.cachedReadGroups());
        assertEquals(1, mReadCache.cachedReadCount());
        assertEquals(1, mReadCache.cachedFragCoordGroups());

        // some more duplicates in the next group
        SAMRecord read7a = createRecord(CHR_1, 150, false);
        SAMRecord read7b = createRecord(CHR_1, 150, false);
        SAMRecord read7c = createRecord(CHR_1, 150, false);

        mReadCache.processRead(read7a);
        mReadCache.processRead(read7b);
        mReadCache.processRead(read7c);

        assertEquals(2, mReadCache.cachedReadGroups());
        assertEquals(4, mReadCache.cachedReadCount());
        assertEquals(2, mReadCache.cachedFragCoordGroups());

        fragCoordReads = mReadCache.evictAll();
        assertNotNull(fragCoordReads);

        assertEquals(1, fragCoordReads.SingleReads.size());
        assertEquals(1, fragCoordReads.DuplicateGroups.size());
        assertEquals(3, fragCoordReads.DuplicateGroups.get(0).readCount());

        assertEquals(0, mReadCache.cachedReadGroups());
        assertEquals(0, mReadCache.cachedReadCount());
        assertEquals(0, mReadCache.cachedFragCoordGroups());

    }

    private static SAMRecord createRecord(final String chromosome, final int readStart, boolean isReversed)
    {
        return createSamRecord(
                READ_ID_GEN.nextId(), chromosome, readStart, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 200, isReversed, false, null, true, TEST_READ_CIGAR);
    }
}
