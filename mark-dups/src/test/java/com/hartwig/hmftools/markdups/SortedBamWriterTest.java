package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_ID;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.markdups.write.SortedBamWriter;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SortedBamWriterTest
{
    @Test
    public void testSortedReadCache()
    {
        SortedBamWriter readCache = new SortedBamWriter(10, 50, null);

        readCache.initialiseStartPosition(CHR_1, 1);
        readCache.setUpperBoundPosition(110);

        readCache.addRecord(createRead(CHR_1, 100));
        readCache.addRecord(createRead(CHR_1, 110));

        assertEquals(0, readCache.written());
        assertEquals(2, readCache.cached());

        // add a record past the buffer position
        readCache.setUpperBoundPosition(170);
        readCache.addRecord(createRead(CHR_1, 170));

        assertEquals(2, readCache.written());
        assertEquals(1, readCache.cached());

        readCache.addRecord(createRead(CHR_1, 180));

        assertEquals(2, readCache.written());
        assertEquals(2, readCache.cached());

        // grow capacity for records within the position buffer
        for(int i = 0; i < 10; ++i)
        {
            readCache.addRecord(createRead(CHR_1, 181));
        }

        assertEquals(2, readCache.written());
        assertEquals(12, readCache.cached());
        // assertEquals(20, readCache.capacity());

        readCache.flush();

        // new chromosome
        readCache.addRecord(createRead(CHR_2, 50));

        readCache.initialiseStartPosition(CHR_2, 1);
        readCache.setUpperBoundPosition(50);

        assertEquals(14, readCache.written());
        assertEquals(1, readCache.cached());
        assertEquals(10, readCache.capacity());
        assertEquals(2, readCache.writeCount());
    }

    private static SAMRecord createRead(final String chromosome, final int position)
    {
        return SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, chromosome, position, "", TEST_READ_CIGAR, CHR_1, 100, false,
                false, null, false, TEST_READ_CIGAR);
    }
}
