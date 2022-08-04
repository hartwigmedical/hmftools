package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.svprep.SvPrepTestUtils.CHR_1;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.createSamRecord;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.readIdStr;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.ExpectedRead;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadRecord;
import com.hartwig.hmftools.svprep.reads.ReadType;

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
        Map<String, List<ExpectedRead>> missedReadsMap = Maps.newHashMap();

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
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap, missedReadsMap);
        assertFalse(read1.written());
        assertEquals(2, getExpectedReadsCount(null));
        assertEquals(1, mCombinedReadGroups.getExpectedReadIds(REGION_2).size());

        spanningGroupsMap.clear();
        missedReadsMap.clear();
        spanningGroupsMap.put(rg2.id(), rg2);
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap, missedReadsMap);
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
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap, missedReadsMap);
        assertFalse(read1.written());
        assertFalse(read2.written());
        assertEquals(2, getExpectedReadsCount(null));

        spanningGroupsMap.clear();
        missedReadsMap.clear();
        spanningGroupsMap.put(rg2.id(), rg2);
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap, missedReadsMap);
        assertTrue(read3.written());
        assertTrue(read4.written());
        assertTrue(missedReadsMap.isEmpty());

        assertEquals(0, getExpectedReadsCount(null));

        // now with the second partition getting the missed read from the first
        mCombinedReadGroups.reset();

        read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));

        rg1 = new ReadGroup(read1);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);

        assertEquals(2, rg1.partitionCount());

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap, missedReadsMap);
        assertFalse(read1.written());
        assertEquals(2, getExpectedReadsCount(null));

        // an empty partition
        spanningGroupsMap.clear();
        missedReadsMap.clear();
        assertEquals(1, mCombinedReadGroups.getExpectedReadIds(REGION_2).size());
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap, missedReadsMap);
        assertEquals(1, missedReadsMap.size());
        ExpectedRead missedRead = missedReadsMap.get(rg1.id()).get(0);
        assertEquals(10800, missedRead.Position);
        assertEquals(false, missedRead.found());

        assertEquals(0, getExpectedReadsCount(null));
    }

    @Test
    public void testCombinedGroupsSupplementary()
    {
        int readId = 0;

        Map<String, ReadGroup> spanningGroupsMap = Maps.newHashMap();
        Map<String, List<ExpectedRead>> missedReadsMap = Maps.newHashMap();

        // 2 reads, one with a supp, spanning different partitions
        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 801, true, false, ""));

        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 801, CHR_1, 800, false, false, "1;10800;-;46S30M;255;0"));

        ReadGroup rg1 = new ReadGroup(read1);
        rg1.addRead(read2);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap, missedReadsMap);
        assertEquals(3, getExpectedReadsCount(null));

        ReadRecord read3 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, true, "1;801;-;46S30M;255;0"));

        rg1 = new ReadGroup(read3);
        rg1.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());
        assertTrue(rg1.conditionalOnRemoteReads());

        spanningGroupsMap.clear();
        missedReadsMap.clear();

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap, missedReadsMap);
        assertEquals(0, getExpectedReadsCount(null));
        assertTrue(rg1.hasRemoteJunctionReads());

        mCombinedReadGroups.reset();
        spanningGroupsMap.clear();
        missedReadsMap.clear();

        // repeat but in the reverse order
        rg1 = new ReadGroup(read3);
        rg1.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap, missedReadsMap);
        assertEquals(3, getExpectedReadsCount(null));

        rg1 = new ReadGroup(read1);
        rg1.addRead(read2);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());

        spanningGroupsMap.clear();
        missedReadsMap.clear();

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap, missedReadsMap);
        assertEquals(0, getExpectedReadsCount(null));
        assertEquals(3, rg1.reads().size());

        // test a group spanning 3 partitions
        mCombinedReadGroups.reset();
        spanningGroupsMap.clear();
        missedReadsMap.clear();

        read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10801, true, false, ""));

        rg1 = new ReadGroup(read1);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap, missedReadsMap);
        assertEquals(2, getExpectedReadsCount(null));

        read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10801, CHR_1, 800, false, false, "1;20800;-;46S30M;255;0"));

        rg1 = new ReadGroup(read2);
        rg1.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(3, rg1.partitionCount());

        spanningGroupsMap.clear();
        missedReadsMap.clear();
        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap, missedReadsMap);
        assertEquals(3, getExpectedReadsCount(null));

        // now the supplementary
        read3 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 20800, CHR_1, 800, false, true, "1;10801;-;46S30M;255;0"));

        rg1 = new ReadGroup(read3);
        rg1.setPartitionCount(REGION_3, PARTITION_SIZE);
        assertEquals(3, rg1.partitionCount());

        spanningGroupsMap.clear();
        missedReadsMap.clear();
        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_3, spanningGroupsMap, missedReadsMap);
        assertEquals(0, getExpectedReadsCount(null));
        assertTrue(rg1.hasRemoteJunctionReads());
    }

    @Test
    public void testCombinedGroupsRemoteCandidates()
    {
        int readId = 0;

        Map<String, ReadGroup> spanningGroupsMap = Maps.newHashMap();
        Map<String, List<ExpectedRead>> missedReadsMap = Maps.newHashMap();

        // a read group split across partitions, both only candidates for support, neither should be written
        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));

        read1.setReadType(ReadType.CANDIDATE_SUPPORT);

        ReadGroup rg1 = new ReadGroup(read1);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());
        assertTrue(rg1.conditionalOnRemoteReads());

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap, missedReadsMap);
        assertEquals(2, getExpectedReadsCount(null));

        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, false, ""));

        read2.setReadType(ReadType.CANDIDATE_SUPPORT);

        ReadGroup rg2 = new ReadGroup(read2);
        rg2.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, rg2.partitionCount());
        assertTrue(rg2.conditionalOnRemoteReads());

        spanningGroupsMap.clear();
        missedReadsMap.clear();

        spanningGroupsMap.put(rg2.id(), rg2);
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap, missedReadsMap);
        assertEquals(0, getExpectedReadsCount(null));
        assertFalse(rg2.hasRemoteJunctionReads()); // won't be written to the BAM for this reason
        assertTrue(missedReadsMap.isEmpty());
        assertEquals(1, rg2.reads().size());

        mCombinedReadGroups.reset();
        spanningGroupsMap.clear();
        missedReadsMap.clear();

        // this time with the second group containing an actual junction fragment
        read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));

        read1.setReadType(ReadType.CANDIDATE_SUPPORT);

        rg1 = new ReadGroup(read1);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());
        assertTrue(rg1.conditionalOnRemoteReads());

        spanningGroupsMap.put(rg1.id(), rg1);
        mCombinedReadGroups.processSpanningReadGroups(REGION_1, spanningGroupsMap, missedReadsMap);
        assertEquals(2, getExpectedReadsCount(null));

        read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, false, ""));

        read2.setReadType(ReadType.JUNCTION);

        rg2 = new ReadGroup(read2);
        rg2.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, rg2.partitionCount());

        spanningGroupsMap.clear();
        missedReadsMap.clear();

        spanningGroupsMap.put(rg2.id(), rg2);
        mCombinedReadGroups.processSpanningReadGroups(REGION_2, spanningGroupsMap, missedReadsMap);
        assertEquals(0, getExpectedReadsCount(null));
        assertEquals(2, rg2.reads().size()); // has picked up the cached candidate support read

    }
}
