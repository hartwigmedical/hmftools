package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.svprep.SvPrepTestUtils.CHR_1;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.REGION_1;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.REGION_2;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.REGION_3;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.createSamRecord;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.readIdStr;
import static com.hartwig.hmftools.svprep.reads.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.JUNCTION;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadRecord;
import com.hartwig.hmftools.svprep.reads.ReadType;

import org.junit.Test;

public class SpanningReadsCacheTest
{
    private final SvConfig mConfig;
    private final SpanningReadCache mSpanningReadCache;

    public SpanningReadsCacheTest()
    {
        mConfig = new SvConfig(PARTITION_SIZE);
        mSpanningReadCache = new SpanningReadCache(mConfig);
    }

    private int getCachedReadsCount(final String readId)
    {
        return mSpanningReadCache.getCachedReadsCount(readId);
        /*
        return mSpanningReadCache.chrPartitionReadGroupsMap().values().stream()
                .filter(x -> readId == null || x.containsKey(readId))
                .mapToInt(x -> x.values().stream().mapToInt(y -> y.size()).sum()).sum();
        */
    }

    private boolean hasJunctionReadId(final String readId)
    {
        return mSpanningReadCache.junctionPartitionReadIdsMap().values().stream()
                .anyMatch(x -> x.contains(readId));
    }

    @Test
    public void testBasicReads()
    {
        int readId = 0;

        Map<String,ReadGroup> spanningGroupsMap = Maps.newHashMap();

        // test 1: 2 junction reads
        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));
        read1.setReadType(JUNCTION);

        ReadGroup readGroup = new ReadGroup(read1);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, readGroup.partitionCount());

        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, false, ""));
        read2.setReadType(JUNCTION);

        readGroup = new ReadGroup(read2);
        readGroup.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertEquals(0, getCachedReadsCount(null));
        assertFalse(hasJunctionReadId(readGroup.id()));
        // assertEquals(1, mSpanningReadCache.getExpectedReadIds(REGION_2).size());

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_2, spanningGroupsMap);

        assertFalse(hasJunctionReadId(readGroup.id()));
        assertEquals(0, getCachedReadsCount(null));

        mSpanningReadCache.reset();


        // test 2: junction then a candidate read
        read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));
        read1.setReadType(JUNCTION);

        readGroup = new ReadGroup(read1);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);

        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);

        read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, false, ""));
        read2.setReadType(CANDIDATE_SUPPORT);

        readGroup = new ReadGroup(read2);
        readGroup.setPartitionCount(REGION_2, PARTITION_SIZE);

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertEquals(0, getCachedReadsCount(null));
        assertFalse(hasJunctionReadId(readGroup.id()));

        mSpanningReadCache.reset();


        // test 3: candidate then junction read
        read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));
        read1.setReadType(CANDIDATE_SUPPORT);

        readGroup = new ReadGroup(read1);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);

        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertEquals(1, getCachedReadsCount(null));
        assertFalse(hasJunctionReadId(readGroup.id()));

        read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, false, ""));
        read2.setReadType(ReadType.SUPPORT);

        readGroup = new ReadGroup(read2);
        readGroup.setPartitionCount(REGION_2, PARTITION_SIZE);

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertEquals(0, getCachedReadsCount(null));
        assertFalse(hasJunctionReadId(readGroup.id()));
        assertTrue(readGroup.reads().contains(read1));

        mSpanningReadCache.reset();


        // test 4: candidate then an empty partition clears out its cache
        read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));
        read1.setReadType(CANDIDATE_SUPPORT);

        readGroup = new ReadGroup(read1);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertEquals(1, getCachedReadsCount(null));

        spanningGroupsMap.clear();
        mSpanningReadCache.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertEquals(0, getCachedReadsCount(null));
    }

    @Test
    public void testCombinedGroupsSupplementary()
    {
        int readId = 0;

        Map<String, ReadGroup> spanningGroupsMap = Maps.newHashMap();

        // test1: 2 reads, one with a supp, spanning different partitions
        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 801, true, false, ""));

        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 801, CHR_1, 800, false, false, "1;10800;-;46S30M;255;0"));
        read2.setReadType(JUNCTION);

        ReadGroup readGroup = new ReadGroup(read1);
        readGroup.addRead(read2);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertTrue(hasJunctionReadId(readGroup.id()));
        assertEquals(0, getCachedReadsCount(null));

        ReadRecord read3 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, true, "1;801;-;46S30M;255;0"));
        read2.setReadType(CANDIDATE_SUPPORT);

        readGroup = new ReadGroup(read3);
        readGroup.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.clear();

        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertEquals(0, getCachedReadsCount(null));
        assertFalse(hasJunctionReadId(readGroup.id()));

        mSpanningReadCache.reset();
        spanningGroupsMap.clear();

        // test 2: a group spanning 3 partitions
        mSpanningReadCache.reset();
        spanningGroupsMap.clear();

        read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10801, true, false, ""));
        read1.setReadType(CANDIDATE_SUPPORT);

        readGroup = new ReadGroup(read1);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertEquals(1, getCachedReadsCount(null));

        read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10801, CHR_1, 800, false, false, "1;20800;-;46S30M;255;0"));
        read2.setReadType(JUNCTION);

        readGroup = new ReadGroup(read2);
        readGroup.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(3, readGroup.partitionCount());

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertEquals(0, getCachedReadsCount(null));
        assertTrue(hasJunctionReadId(readGroup.id()));  // keeps the junction id for expected reads
        assertTrue(readGroup.reads().contains(read1)); // has picked up the candidate

        // now the supplementary
        read3 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 20800, CHR_1, 800, false, true, "1;10801;-;46S30M;255;0"));
        read3.setReadType(CANDIDATE_SUPPORT);

        readGroup = new ReadGroup(read3);
        readGroup.setPartitionCount(REGION_3, PARTITION_SIZE);
        assertEquals(3, readGroup.partitionCount());

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_3, spanningGroupsMap);
        assertEquals(0, getCachedReadsCount(null));
        assertFalse(hasJunctionReadId(readGroup.id()));


        // test 3: a supplementary and a candidate, followed by the junction
        mSpanningReadCache.reset();
        spanningGroupsMap.clear();

        read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10801, false, false, ""));
        read1.setReadType(CANDIDATE_SUPPORT);

        read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 801, CHR_1, 800, true, true, "1;10801;-;46S30M;255;0"));
        read2.setReadType(JUNCTION);

        readGroup = new ReadGroup(read1);
        readGroup.addRead(read2);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);

        // group treated as a junction
        assertEquals(0, getCachedReadsCount(null));
        assertTrue(hasJunctionReadId(readGroup.id())); // doesn't know if the following read will be a junction

        // now the junction
        read3 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10801, CHR_1, 800, true, false, "1;801;-;46S30M;255;0"));
        read3.setReadType(JUNCTION);

        readGroup = new ReadGroup(read3);
        readGroup.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertEquals(0, getCachedReadsCount(null));
        assertFalse(hasJunctionReadId(readGroup.id()));
        assertEquals(1, readGroup.reads().size());
    }

    /*
    @Test
    public void testSupplementaryLast()
    {
        int readId = 0;

        Map<String, ReadGroup> spanningGroupsMap = Maps.newHashMap();

        // test3: the 2 candidates first then the junction
        mSpanningReadCache.reset();
        spanningGroupsMap.clear();

        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10801, true, false, "1;20800;-;46S30M;255;0"));
        read1.setReadType(CANDIDATE_SUPPORT);

        ReadGroup readGroup = new ReadGroup(read1);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);

        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertEquals(1, getCachedReadsCount(null));

        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10801, CHR_1, 800, false, false, ""));
        read2.setReadType(CANDIDATE_SUPPORT);

        readGroup = new ReadGroup(read2);
        readGroup.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        // assertEquals(1, getCachedReadsCount(null));
        assertEquals(2, getCachedReadsCount(null)); // doesn't know about the supp read to come, so isn't stored

        // now the supplementary
        ReadRecord read3 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 20800, CHR_1, 10801, true, true, "1;800;-;46S30M;255;0"));
        read3.setReadType(JUNCTION);

        readGroup = new ReadGroup(read3);
        readGroup.setPartitionCount(REGION_3, PARTITION_SIZE);

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_3, spanningGroupsMap);
        assertEquals(0, getCachedReadsCount(null));
        assertFalse(hasJunctionReadId(readGroup.id()));
        assertTrue(readGroup.reads().contains(read1)); // has picked up 2 candidates
        assertTrue(readGroup.reads().contains(read2));
    }
    */

    @Test
    public void testCandidateOnlyGroups()
    {
        int readId = 0;

        Map<String, ReadGroup> spanningGroupsMap = Maps.newHashMap();

        // test 1: a read group split across partitions, both only candidates for support, neither should be written
        ReadRecord read1 = ReadRecord.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));

        read1.setReadType(CANDIDATE_SUPPORT);

        ReadGroup rg1 = new ReadGroup(read1);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());
        assertTrue(rg1.conditionalOnRemoteReads());

        spanningGroupsMap.put(rg1.id(), rg1);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertEquals(1, getCachedReadsCount(null));

        ReadRecord read2 = ReadRecord.from(createSamRecord(
                readIdStr(readId), CHR_1, 10800, CHR_1, 800, false, false, ""));

        read2.setReadType(CANDIDATE_SUPPORT);

        ReadGroup rg2 = new ReadGroup(read2);
        rg2.setPartitionCount(REGION_2, PARTITION_SIZE);
        assertEquals(2, rg2.partitionCount());
        assertTrue(rg2.conditionalOnRemoteReads());

        spanningGroupsMap.clear();

        spanningGroupsMap.put(rg2.id(), rg2);
        mSpanningReadCache.processSpanningReadGroups(REGION_2, spanningGroupsMap);
        assertEquals(0, getCachedReadsCount(null));
        assertFalse(rg2.hasRemoteJunctionReads());
        assertEquals(1, rg2.reads().size());

        mSpanningReadCache.reset();
        spanningGroupsMap.clear();
    }
}
