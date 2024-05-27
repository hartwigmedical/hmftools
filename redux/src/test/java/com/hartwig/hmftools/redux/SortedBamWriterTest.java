package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_ID;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.write.SortedBamConfig;
import com.hartwig.hmftools.redux.write.SortedBamWriter;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SortedBamWriterTest
{
    @Test
    public void testSortedReadCache()
    {
        SortedBamConfig config = new SortedBamConfig(20, 5, 1);

        SortedBamWriter readCache = new SortedBamWriter(config, null);

        readCache.initialiseStartPosition(CHR_1, 50);
        readCache.setUpperWritablePosition(100);
        readCache.setUpperBoundPosition(110);

        // first test boundaries
        assertTrue(readCache.canWriteRecord(createRead(CHR_1, 105)));
        assertFalse(readCache.canWriteRecord(createRead(CHR_2, 105))); // wrong chromosome
        assertFalse(readCache.canWriteRecord(createRead(CHR_1, 115))); // too high

        readCache.addRecord(createRead(CHR_1, 100));
        readCache.addRecord(createRead(CHR_1, 110));

        assertEquals(0, readCache.written());
        assertEquals(2, readCache.cached());

        // add a record past the buffer position
        readCache.setUpperBoundPosition(170);
        readCache.setUpperWritablePosition(120);
        readCache.addRecord(createRead(CHR_1, 170));

        assertEquals(2, readCache.written());
        assertEquals(1, readCache.cached());

        assertFalse(readCache.canWriteRecord(createRead(CHR_1, 40))); // too low

        readCache.addRecord(createRead(CHR_1, 180));

        assertEquals(2, readCache.written());
        assertEquals(2, readCache.cached());

        // grow capacity for records within the position buffer
        readCache.setUpperBoundPosition(200);

        for(int i = 0; i < 10; ++i)
        {
            readCache.addRecord(createRead(CHR_1, 181));
        }

        assertEquals(2, readCache.written());
        assertEquals(12, readCache.cached());

        readCache.setUpperWritablePosition(190);

        readCache.addRecord(createRead(CHR_1, 200));

        assertEquals(14, readCache.written());
        assertEquals(1, readCache.cached());

        readCache.flush();

        assertEquals(15, readCache.written());
        assertEquals(0, readCache.cached());

        // new chromosome
        readCache.initialiseStartPosition(CHR_2, 1);
        readCache.setUpperBoundPosition(50);

        readCache.addRecord(createRead(CHR_2, 50));

        assertEquals(15, readCache.written());
        assertEquals(1, readCache.cached());
        assertEquals(3, readCache.writeCount());
    }

    private static SAMRecord createRead(final String chromosome, final int position)
    {
        return SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, chromosome, position, "", TEST_READ_CIGAR, CHR_1, 100, false,
                false, null, false, TEST_READ_CIGAR);
    }
}
