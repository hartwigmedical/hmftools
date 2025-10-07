package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.redux.common.DuplicatesConfig;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UnpairedReadDuplicateTest
{
    public UnpairedReadDuplicateTest()
    {
        ReduxConfig.SEQUENCING_TYPE = SequencingType.SBX;
    }

    @After
    public void resetSequencingType()
    {
        ReduxConfig.SEQUENCING_TYPE = SequencingType.ILLUMINA;
    }

    @Test
    public void testSbxDuplicateGroupCollapsing()
    {
        int maxDuplicateDistance = 2;
        DuplicatesConfig groupCollapseConfig = new DuplicatesConfig(maxDuplicateDistance);
        ReadCache readCache = new ReadCache(100, 100, false, groupCollapseConfig);

        // no collapsing
        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150 + maxDuplicateDistance + 1, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100 - maxDuplicateDistance - 1, 150, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150 - maxDuplicateDistance - 1, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100 + maxDuplicateDistance + 1, 150, false));

        FragmentCoordReads fragmentCoordsReads = readCache.evictAll();

        assertEquals(0, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(5, fragmentCoordsReads.SingleReads.size());
        assertEquals(5, fragmentCoordsReads.totalReadCount());

        // no collapsing of forward and reverse reads
        readCache.clear();

        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 50, 100, true));

        fragmentCoordsReads = readCache.evictAll();

        assertEquals(0, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(2, fragmentCoordsReads.SingleReads.size());
        assertEquals(2, fragmentCoordsReads.totalReadCount());

        // simple collapsing
        readCache.clear();

        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150 + maxDuplicateDistance, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100 - maxDuplicateDistance / 2, 150 + maxDuplicateDistance / 2, false));

        fragmentCoordsReads = readCache.evictAll();

        assertEquals(1, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(0, fragmentCoordsReads.SingleReads.size());
        assertEquals(3, fragmentCoordsReads.totalReadCount());

        // no multi-coord duplicate groups
        readCache.clear();

        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150, false));

        fragmentCoordsReads = readCache.evictAll();

        assertEquals(1, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(0, fragmentCoordsReads.SingleReads.size());
        assertEquals(2, fragmentCoordsReads.totalReadCount());

        // no chain collapsing
        readCache.clear();

        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150 + maxDuplicateDistance, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100 - maxDuplicateDistance, 150 + maxDuplicateDistance, false));
        readCache.processRead(createUnpairedRecord(CHR_1,
                100 - maxDuplicateDistance / 2, 150 + maxDuplicateDistance + maxDuplicateDistance / 2, false));

        fragmentCoordsReads = readCache.evictAll();

        assertEquals(2, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(0, fragmentCoordsReads.SingleReads.size());
        assertEquals(4, fragmentCoordsReads.totalReadCount());

        // a more complex scenario
        readCache.clear();

        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150 + maxDuplicateDistance, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100, 150 - maxDuplicateDistance, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100 + maxDuplicateDistance, 150, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 100 - maxDuplicateDistance, 150, false));

        readCache.processRead(createUnpairedRecord(CHR_1, 100, 200, true));
        readCache.processRead(createUnpairedRecord(CHR_1, 100, 200, true));

        readCache.processRead(createUnpairedRecord(CHR_1, 100, 200, false));

        readCache.processRead(createUnpairedRecord(CHR_1, 200, 300, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 200 + maxDuplicateDistance / 2, 300 + maxDuplicateDistance / 2, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 200 + maxDuplicateDistance / 2, 300 - maxDuplicateDistance / 2, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 200 - maxDuplicateDistance / 2, 300 + maxDuplicateDistance / 2, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 200 - maxDuplicateDistance / 2, 300 - maxDuplicateDistance / 2, false));

        fragmentCoordsReads = readCache.evictAll();

        assertEquals(4, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(4, fragmentCoordsReads.SingleReads.size());
        assertEquals(13, fragmentCoordsReads.totalReadCount());
    }

    @Test
    public void testSbxDuplicateGroupCollapsingLargestGroupFirst()
    {
        int maxDuplicateDistnace = 1;
        DuplicatesConfig groupCollapseConfig = new DuplicatesConfig(maxDuplicateDistnace);
        ReadCache readCache = new ReadCache(100, 100, false, groupCollapseConfig);

        SAMRecord read1 = createUnpairedRecord(CHR_1, 100, 150, false);

        SAMRecord read2 = createUnpairedRecord(CHR_1, 100, 151, false);
        SAMRecord read3 = createUnpairedRecord(CHR_1, 100, 151, false);

        SAMRecord read4 = createUnpairedRecord(CHR_1, 100, 152, false);
        SAMRecord read5 = createUnpairedRecord(CHR_1, 100, 152, false);
        SAMRecord read6 = createUnpairedRecord(CHR_1, 100, 152, false);

        readCache.processRead(read1);
        readCache.processRead(read2);
        readCache.processRead(read3);
        readCache.processRead(read4);
        readCache.processRead(read5);
        readCache.processRead(read6);

        FragmentCoordReads fragmentCoordsReads = readCache.evictAll();

        assertEquals(1, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(5, fragmentCoordsReads.DuplicateGroups.get(0).readCount());
        assertEquals(1, fragmentCoordsReads.SingleReads.size());
        assertEquals(6, fragmentCoordsReads.totalReadCount());
    }

    private static SAMRecord createUnpairedRecord(final String chromosome, final int readStart, int readEnd, boolean isReversed)
    {
        return TestUtils.createUnpairedRecord(READ_ID_GEN.nextId(), chromosome, readStart, readEnd, isReversed);
    }
}
