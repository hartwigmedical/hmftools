package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.TestUtils.readIdStr;
import static com.hartwig.hmftools.esvee.prep.TestUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.esvee.prep.TestUtils.REGION_1;
import static com.hartwig.hmftools.esvee.prep.TestUtils.REGION_2;
import static com.hartwig.hmftools.esvee.prep.TestUtils.REGION_3;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.JUNCTION;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadType;

import org.junit.Test;

public class SpanningReadsCacheTest
{
    private final PrepConfig mConfig;
    private final SpanningReadCache mSpanningReadCache;

    public SpanningReadsCacheTest()
    {
        mConfig = new PrepConfig(PARTITION_SIZE);
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
        PrepRead read1 = PrepRead.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));
        read1.setReadType(JUNCTION);

        ReadGroup readGroup = new ReadGroup(read1);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, readGroup.partitionCount());

        PrepRead read2 = PrepRead.from(createSamRecord(
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
        read1 = PrepRead.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));
        read1.setReadType(JUNCTION);

        readGroup = new ReadGroup(read1);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);

        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.clear();
        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);

        read2 = PrepRead.from(createSamRecord(
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
        read1 = PrepRead.from(createSamRecord(
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

        read2 = PrepRead.from(createSamRecord(
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
        read1 = PrepRead.from(createSamRecord(
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
        PrepRead read1 = PrepRead.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 801, true, false, ""));

        PrepRead read2 = PrepRead.from(createSamRecord(
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

        PrepRead read3 = PrepRead.from(createSamRecord(
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

        read1 = PrepRead.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10801, true, false, ""));
        read1.setReadType(CANDIDATE_SUPPORT);

        readGroup = new ReadGroup(read1);
        readGroup.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, readGroup.partitionCount());

        spanningGroupsMap.put(readGroup.id(), readGroup);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertEquals(1, getCachedReadsCount(null));

        read2 = PrepRead.from(createSamRecord(
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
        read3 = PrepRead.from(createSamRecord(
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

        read1 = PrepRead.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10801, false, false, ""));
        read1.setReadType(CANDIDATE_SUPPORT);

        read2 = PrepRead.from(createSamRecord(
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
        read3 = PrepRead.from(createSamRecord(
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


        // test 4: a supplementary and a candidate, followed by the junction
        mSpanningReadCache.reset();
        spanningGroupsMap.clear();

        read1 = PrepRead.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10801, false, false, ""));
        read1.setReadType(CANDIDATE_SUPPORT);

        read2 = PrepRead.from(createSamRecord(
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
        read3 = PrepRead.from(createSamRecord(
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

    @Test
    public void testCandidateOnlyGroups()
    {
        int readId = 0;

        Map<String, ReadGroup> spanningGroupsMap = Maps.newHashMap();

        // test 1: a read group split across partitions, both only candidates for support, neither should be written
        PrepRead read1 = PrepRead.from(createSamRecord(
                readIdStr(++readId), CHR_1, 800, CHR_1, 10800, true, false, ""));

        read1.setReadType(CANDIDATE_SUPPORT);

        ReadGroup rg1 = new ReadGroup(read1);
        rg1.setPartitionCount(REGION_1, PARTITION_SIZE);
        assertEquals(2, rg1.partitionCount());
        assertTrue(rg1.conditionalOnRemoteReads());

        spanningGroupsMap.put(rg1.id(), rg1);
        mSpanningReadCache.processSpanningReadGroups(REGION_1, spanningGroupsMap);
        assertEquals(1, getCachedReadsCount(null));

        PrepRead read2 = PrepRead.from(createSamRecord(
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
