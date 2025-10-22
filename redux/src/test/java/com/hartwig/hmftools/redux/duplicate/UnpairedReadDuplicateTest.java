package com.hartwig.hmftools.redux.duplicate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.duplicate.ReadCache.DEFAULT_POP_DISTANCE_CHECK;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.TestUtils;

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
        int maxDupDistance = 2;
        DuplicatesConfig groupCollapseConfig = new DuplicatesConfig(maxDupDistance);
        ReadCache readCache = new ReadCache(
                100, 100, false, groupCollapseConfig, DEFAULT_POP_DISTANCE_CHECK, 0);

        int readStart = 100;
        int readEnd = 150;

        // no collapsing
        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd + maxDupDistance + 1, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart - maxDupDistance - 1, readEnd, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd - maxDupDistance - 1, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart + maxDupDistance + 1, readEnd, false));

        FragmentCoordReads fragmentCoordsReads = readCache.evictAll();

        assertEquals(0, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(5, fragmentCoordsReads.SingleReads.size());
        assertEquals(5, fragmentCoordsReads.totalReadCount());

        // no collapsing of forward and reverse reads
        readCache.clear();

        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 50, 100, true));

        fragmentCoordsReads = readCache.evictAll();

        assertEquals(0, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(2, fragmentCoordsReads.SingleReads.size());
        assertEquals(2, fragmentCoordsReads.totalReadCount());

        // simple collapsing
        readCache.clear();

        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd + maxDupDistance, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart - maxDupDistance / 2, readEnd + maxDupDistance / 2, false));

        fragmentCoordsReads = readCache.evictAll();

        assertEquals(1, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(0, fragmentCoordsReads.SingleReads.size());
        assertEquals(3, fragmentCoordsReads.totalReadCount());

        // no multi-coord duplicate groups
        readCache.clear();

        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd, false));

        fragmentCoordsReads = readCache.evictAll();

        assertEquals(1, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(0, fragmentCoordsReads.SingleReads.size());
        assertEquals(2, fragmentCoordsReads.totalReadCount());

        // no chain collapsing
        readCache.clear();

        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd + maxDupDistance, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart - maxDupDistance, readEnd + maxDupDistance, false));
        readCache.processRead(createUnpairedRecord(CHR_1,
                readStart - maxDupDistance / 2, readEnd + maxDupDistance + maxDupDistance / 2, false));

        fragmentCoordsReads = readCache.evictAll();

        assertEquals(1, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(4, fragmentCoordsReads.DuplicateGroups.get(0).totalReadCount());
        assertEquals(0, fragmentCoordsReads.SingleReads.size());
        assertEquals(4, fragmentCoordsReads.totalReadCount());

        // a more complex scenario
        readCache.clear();

        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd + maxDupDistance, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart, readEnd - maxDupDistance, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart + maxDupDistance, readEnd, false));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart - maxDupDistance, readEnd, false));

        readCache.processRead(createUnpairedRecord(CHR_1, readStart, 200, true));
        readCache.processRead(createUnpairedRecord(CHR_1, readStart, 200, true));

        readCache.processRead(createUnpairedRecord(CHR_1, readStart, 200, false));

        readCache.processRead(createUnpairedRecord(CHR_1, 200, 300, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 200 + maxDupDistance / 2, 300 + maxDupDistance / 2, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 200 + maxDupDistance / 2, 300 - maxDupDistance / 2, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 200 - maxDupDistance / 2, 300 + maxDupDistance / 2, false));
        readCache.processRead(createUnpairedRecord(CHR_1, 200 - maxDupDistance / 2, 300 - maxDupDistance / 2, false));

        fragmentCoordsReads = readCache.evictAll();

        assertEquals(3, fragmentCoordsReads.DuplicateGroups.size());
        DuplicateGroup duplicateGroup = fragmentCoordsReads.DuplicateGroups.get(0);
        assertEquals(98, duplicateGroup.reads().stream().mapToInt(x -> x.getAlignmentStart()).min().orElse(0));
        assertEquals(102, duplicateGroup.reads().stream().mapToInt(x -> x.getAlignmentStart()).max().orElse(0));
        assertEquals(148, duplicateGroup.reads().stream().mapToInt(x -> x.getAlignmentEnd()).min().orElse(0));
        assertEquals(152, duplicateGroup.reads().stream().mapToInt(x -> x.getAlignmentEnd()).max().orElse(0));

        duplicateGroup = fragmentCoordsReads.DuplicateGroups.get(2);
        assertEquals(199, duplicateGroup.reads().stream().mapToInt(x -> x.getAlignmentStart()).min().orElse(0));
        assertEquals(201, duplicateGroup.reads().stream().mapToInt(x -> x.getAlignmentStart()).max().orElse(0));
        assertEquals(299, duplicateGroup.reads().stream().mapToInt(x -> x.getAlignmentEnd()).min().orElse(0));
        assertEquals(301, duplicateGroup.reads().stream().mapToInt(x -> x.getAlignmentEnd()).max().orElse(0));

        assertEquals(0, fragmentCoordsReads.SingleReads.size());
        assertEquals(13, fragmentCoordsReads.totalReadCount());
    }

    @Test
    public void testSbxDuplicateGroupCollapsingLargestGroupFirst()
    {
        int maxDupDistance = 1;
        DuplicatesConfig groupCollapseConfig = new DuplicatesConfig(maxDupDistance);
        ReadCache readCache = new ReadCache(
                100, 100, false, groupCollapseConfig, DEFAULT_POP_DISTANCE_CHECK, 0);

        int readStart = 100;
        int readEnd = 250;

        // too far from the largest group and does not do chain mering
        SAMRecord read1 = createUnpairedRecord(CHR_1, readStart, readEnd, false);

        SAMRecord read2 = createUnpairedRecord(CHR_1, readStart, readEnd + 1, false);
        SAMRecord read3 = createUnpairedRecord(CHR_1, readStart, readEnd + 1, false);

        SAMRecord read4 = createUnpairedRecord(CHR_1, readStart, readEnd + 2, false);
        SAMRecord read5 = createUnpairedRecord(CHR_1, readStart, readEnd + 2, false);
        SAMRecord read6 = createUnpairedRecord(CHR_1, readStart, readEnd + 2, false);

        readCache.processRead(read1);
        readCache.processRead(read2);
        readCache.processRead(read3);
        readCache.processRead(read4);
        readCache.processRead(read5);
        readCache.processRead(read6);

        FragmentCoordReads fragmentCoordsReads = readCache.evictAll();

        assertEquals(1, fragmentCoordsReads.DuplicateGroups.size());
        assertEquals(6, fragmentCoordsReads.DuplicateGroups.get(0).totalReadCount());
        assertEquals(0, fragmentCoordsReads.SingleReads.size());
        assertEquals(6, fragmentCoordsReads.totalReadCount());
    }

    private static SAMRecord createUnpairedRecord(final String chromosome, final int readStart, int readEnd, boolean isReversed)
    {
        return TestUtils.createUnpairedRecord(READ_ID_GEN.nextId(), chromosome, readStart, readEnd, isReversed);
    }
}
