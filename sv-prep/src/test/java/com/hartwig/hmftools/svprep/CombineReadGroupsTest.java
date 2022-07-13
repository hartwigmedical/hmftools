package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.svprep.SvPrepTestUtils.CHR_1;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.DEFAULT_BASE_QUAL;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.DEFAULT_MAP_QUAL;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.buildFlags;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.createSamRecord;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.readIdStr;
import static com.hartwig.hmftools.svprep.reads.ReadType.JUNCTION;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadRecord;

import org.junit.Test;

public class CombineReadGroupsTest
{
    private static final int PARTITION_SIZE = 10000;
    private static final ChrBaseRegion REGION_1 = new ChrBaseRegion(CHR_1, 1, PARTITION_SIZE - 1);
    private static final ChrBaseRegion REGION_2 = new ChrBaseRegion(CHR_1, REGION_1.end() + 1, REGION_1.end() + PARTITION_SIZE);
    private static final ChrBaseRegion REGION_3 = new ChrBaseRegion(CHR_1, REGION_2.end() + 1, REGION_2.end() + PARTITION_SIZE);

    private final SvConfig mConfig;
    private final CombinedReadGroups mCombinedReadGroups;

    public CombineReadGroupsTest()
    {
        mConfig = new SvConfig(PARTITION_SIZE);
        mCombinedReadGroups = new CombinedReadGroups(mConfig);
    }

    private int getExpectedReadsCount(final String readId)
    {
        return mCombinedReadGroups.chrPartitionReadGroupsMap().values().stream()
                .filter(x -> readId == null || x.containsKey(readId))
                .mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();
    }

    @Test
    public void testCombinedGroupsBasic()
    {
        int readId = 0;

        Map<String, ReadGroup> spanningGroupsMap = Maps.newHashMap();

        // 2 reads, no supp, spanning different partitions
        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));

        ReadGroup rg1 = new ReadGroup(read1);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());

        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, false, ""));

        ReadGroup rg2 = new ReadGroup(read2);
        rg2.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, rg2.partitionCount());

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertFalse(read1.written());
        assertEquals(2, getExpectedReadsCount(null));

        spanningGroupsMap.clear();
        spanningGroupsMap.put(rg2.id(), rg2);
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertFalse(read2.written());

        assertEquals(0, getExpectedReadsCount(null));

        // same again but with each partition getting the other's read
        mCombinedReadGroups.reset();

        read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));

        read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, false, ""));

        rg1 = new ReadGroup(read1);
        rg1.addRead(read2);

        ReadRecord read3 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 800, CHR_1, 10800, true, false, ""));

        ReadRecord read4 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, false, ""));

        rg2 = new ReadGroup(read4);
        rg2.addRead(read3);

        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        rg2.setPartitionCount(REGION_2, PARTITION_SIZE);

        assertEquals(2, rg1.partitionCount());
        assertEquals(2, rg2.partitionCount());

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertFalse(read1.written());
        assertFalse(read2.written());
        assertEquals(2, getExpectedReadsCount(null));

        spanningGroupsMap.clear();
        spanningGroupsMap.put(rg2.id(), rg2);
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertTrue(read3.written());
        assertTrue(read4.written());

        assertEquals(0, getExpectedReadsCount(null));
    }

    @Test
    public void testCombinedGroupsSupplementary()
    {
        int readId = 0;

        Map<String, ReadGroup> spanningGroupsMap = Maps.newHashMap();

        // 2 reads, no supp, spanning different partitions
        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 801, true, false, ""));

        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 801, CHR_1, 800, false, false, "1;10800;-;46S30M;255;0"));

        ReadGroup rg1 = new ReadGroup(read1);
        rg1.addRead(read2);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(1, rg1.partitionCount()); // supplementaries aren't counted towards remote partitions

        /*
        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertEquals(3, getExpectedReadsCount(null));

        ReadRecord read3 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, true, "1;801;-;46S30M;255;0"));

        ReadGroup rg2 = new ReadGroup(read3);
        rg2.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, rg2.partitionCount());

        spanningGroupsMap.clear();
        spanningGroupsMap.put(rg2.id(), rg2);
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertFalse(read1.written());
        assertFalse(read2.written());
        assertFalse(read3.written());

        assertEquals(0, getExpectedReadsCount(null));
        */

        // test a group spanning 3 partitions (though the supp is ignored)
        mCombinedReadGroups.reset();
        spanningGroupsMap.clear();

        read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10801, true, false, ""));

        read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10801, CHR_1, 800, false, false, "1;20800;-;46S30M;255;0"));

        rg1 = new ReadGroup(read1);
        rg1.addRead(read2);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertEquals(2, getExpectedReadsCount(null));

        /*
        ReadRecord read3 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 20800, CHR_1, 800, false, true, "1;10801;-;46S30M;255;0"));

        rg2 = new ReadGroup(read3);
        rg2.setPartitionCount(REGION_3, PARTITION_SIZE);
        assertEquals(3, rg2.partitionCount());

        spanningGroupsMap.clear();
        spanningGroupsMap.put(rg2.id(), rg2);
        mCombinedReadGroups.processSpanningReadGroups(REGION_3, spanningGroupsMap);
        assertFalse(read1.written());
        assertFalse(read2.written());
        assertFalse(read3.written());

        assertEquals(1, getExpectedReadsCount(null));
        */

        // send partition misses the second read and so doesn't get its mate from the first partition either
        spanningGroupsMap.clear();
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertEquals(0, getExpectedReadsCount(null));
    }
}
