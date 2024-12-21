package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadCacheTest
{
    @Test
    public void testBasics()
    {
        ReadCache readCache = new ReadCache(100, 100, false);

        SAMRecord read1 = createRecord(CHR_1, 50, false);

        readCache.processRead(read1);

        assertEquals(1, readCache.cachedReadGroups());
        assertEquals(1, readCache.cachedFragCoordGroups());
        assertEquals(1, readCache.cachedReadCount());

        FragmentCoordReads fragCoordReads = readCache.popReads();
        assertNull(fragCoordReads);

        SAMRecord read2a = createRecord(CHR_1, 80, false);
        SAMRecord read2b = createRecord(CHR_1, 80, false);

        readCache.processRead(read2a);
        readCache.processRead(read2b);

        SAMRecord read3 = createRecord(CHR_1, 120, false);

        readCache.processRead(read3);

        assertEquals(2, readCache.cachedReadGroups());
        assertEquals(4, readCache.cachedReadCount());
        assertEquals(3, readCache.cachedFragCoordGroups());

        fragCoordReads = readCache.popReads();
        assertNull(fragCoordReads);

        SAMRecord read4 = createRecord(CHR_1, 199, false);

        readCache.processRead(read4);

        assertEquals(2, readCache.cachedReadGroups());

        // the next read triggers the first group to be processed
        SAMRecord read5 = createRecord(CHR_1, 201, false);

        readCache.processRead(read5);

        fragCoordReads = readCache.popReads();
        assertNotNull(fragCoordReads);

        assertEquals(1, fragCoordReads.SingleReads.size());
        assertEquals(1, fragCoordReads.DuplicateGroups.size());

        assertEquals(2, readCache.cachedReadGroups());
        assertEquals(3, readCache.cachedReadCount());
        assertEquals(3, readCache.cachedFragCoordGroups());

        // next read is on a new chromosome so clears all previous
        SAMRecord read6 = createRecord(CHR_2, 50, false);

        readCache.processRead(read6);

        fragCoordReads = readCache.popReads();
        assertNotNull(fragCoordReads);

        assertEquals(3, fragCoordReads.SingleReads.size());
        assertEquals(0, fragCoordReads.DuplicateGroups.size());

        assertEquals(1, readCache.cachedReadGroups());
        assertEquals(1, readCache.cachedReadCount());
        assertEquals(1, readCache.cachedFragCoordGroups());

        // some more duplicates in the next group
        SAMRecord read7a = createRecord(CHR_1, 150, false);
        SAMRecord read7b = createRecord(CHR_1, 150, false);
        SAMRecord read7c = createRecord(CHR_1, 150, false);

        readCache.processRead(read7a);
        readCache.processRead(read7b);
        readCache.processRead(read7c);

        assertEquals(2, readCache.cachedReadGroups());
        assertEquals(4, readCache.cachedReadCount());
        assertEquals(2, readCache.cachedFragCoordGroups());

        fragCoordReads = readCache.evictAll();
        assertNotNull(fragCoordReads);

        assertEquals(1, fragCoordReads.SingleReads.size());
        assertEquals(1, fragCoordReads.DuplicateGroups.size());
        assertEquals(3, fragCoordReads.DuplicateGroups.get(0).readCount());

        assertEquals(0, readCache.cachedReadGroups());
        assertEquals(0, readCache.cachedReadCount());
        assertEquals(0, readCache.cachedFragCoordGroups());
    }

    @Test
    public void testMixedPrimarySupplementaries()
    {
        ReadCache readCache = new ReadCache(100, 100, false);

        SupplementaryReadData suppData1 = new SupplementaryReadData(CHR_2, 100, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

        SAMRecord read1 = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 200, false, false, suppData1, true, TEST_READ_CIGAR);

        SupplementaryReadData suppData2 = new SupplementaryReadData(CHR_2, 100, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

        SAMRecord read2 = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 200, false, true, suppData2, true, TEST_READ_CIGAR);

        readCache.processRead(read1);
        readCache.processRead(read2);

        assertEquals(1, readCache.cachedReadGroups());
        assertEquals(2, readCache.cachedReadCount());
        assertEquals(2, readCache.cachedFragCoordGroups());
    }

    private static SAMRecord createRecord(final String chromosome, final int readStart, boolean isReversed)
    {
        return createSamRecord(
                READ_ID_GEN.nextId(), chromosome, readStart, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 200, isReversed, false, null, true, TEST_READ_CIGAR);
    }
}
