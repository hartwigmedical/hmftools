package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.flipFirstInPair;
import static com.hartwig.hmftools.redux.common.DuplicateGroupCollapser.SINGLE_END_JITTER_COLLAPSE_DISTANCE;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;
import com.hartwig.hmftools.redux.common.ReadInfo;

import org.junit.jupiter.api.Test;

import htsjdk.samtools.SAMRecord;

public class JitterReadCacheTest
{
    private static final String TEST_READ_BASES = "A".repeat(100);
    private static final String TEST_CIGAR = "100M";

    @Test
    public void testHandlesReadCachePossibleHoldingBackOfReadsForPerformanceReasons()
    {
        JitterReadCache readCache =  new JitterReadCache(
                new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, true, ILLUMINA));
        String readBases = "A".repeat(143);
        String cigar = "143M";

        // read1 and read2 need to be popped off together
        SAMRecord read1 = createSamRecord("READ_001", CHR_1, 850, readBases, cigar, CHR_1, 1047, false, false, null, true, cigar);
        SAMRecord read2 = createSamRecord("READ_002", CHR_1, 851, readBases, cigar, CHR_1, 1047, false, false, null, true, cigar);
        SAMRecord read3 = createSamRecord("READ_003", CHR_1, 1000, readBases, cigar, CHR_1, 919, true, false, null, false, cigar);
        SAMRecord read4 = createSamRecord("READ_004", CHR_1, 1098, readBases, cigar, CHR_1, 914, true, false, null, false, cigar);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3, read4);

        final Function<FragmentCoordReads, Stream<SAMRecord>> fragmentCoordReadsToReads =
                (final FragmentCoordReads fragmentCoordReads) -> Stream.concat(
                        fragmentCoordReads.DuplicateGroups.stream().flatMap(x -> x.reads().stream()),
                        fragmentCoordReads.SingleReads.stream().map(ReadInfo::read)
                );

        Set<String> readNamesOfInterest = Sets.newHashSet("READ_001", "READ_002");
        int totalPopsOfInterest = 0;
        for(SAMRecord read : reads)
        {
            readCache.processRead(read);
            FragmentCoordReads fragmentCoordReads = readCache.popReads();
            if(fragmentCoordReads == null)
                continue;

            int popsOfInterest = (int) fragmentCoordReadsToReads.apply(fragmentCoordReads)
                    .filter(x -> readNamesOfInterest.contains(x.getReadName()))
                    .count();
            totalPopsOfInterest += popsOfInterest;

            assertTrue(popsOfInterest == 0 || popsOfInterest == 2);
        }

        FragmentCoordReads fragmentCoordReads = readCache.evictAll();
        if(fragmentCoordReads != null)
        {
            int popsOfInterest = (int) fragmentCoordReadsToReads.apply(fragmentCoordReads)
                    .filter(x -> readNamesOfInterest.contains(x.getReadName()))
                    .count();
            totalPopsOfInterest += popsOfInterest;

            assertTrue(popsOfInterest == 0 || popsOfInterest == 2);
        }

        assertEquals(2, totalPopsOfInterest);
    }

    @Test
    public void testReadFlows()
    {
        JitterReadCache readCache = new JitterReadCache(
                new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, false, ILLUMINA));

        assertEquals(-1, readCache.minCachedReadStart());
        assertEquals(0, readCache.cachedReadCount());
        assertEquals(0, readCache.cachedFragCoordGroups());
        assertEquals(0, readCache.currentReadMinPosition());

        final BiFunction<Integer, Integer, SAMRecord> createRead = (final Integer id, final Integer pos) ->
                createSamRecord(format("READ_%03d", id), CHR_1, pos, TEST_READ_BASES, TEST_CIGAR, CHR_2, 100, false, false, null, true, TEST_CIGAR);

        readCache.processRead(createRead.apply(1, 100));
        assertEquals(100, readCache.currentReadMinPosition());
        assertEquals(100, readCache.minCachedReadStart());
        assertEquals(1, readCache.cachedReadCount());
        assertEquals(1, readCache.cachedFragCoordGroups());
        Multiset<String> poppedReads = collaseFragmentCoordReads(readCache.popReads());
        assertTrue(poppedReads.isEmpty());
        assertTrue(readCache.auxCacheReadNames().isEmpty());
        assertEquals(100, readCache.minCachedReadStart());
        assertEquals(1, readCache.cachedReadCount());
        assertEquals(1, readCache.cachedFragCoordGroups());

        readCache.processRead(createRead.apply(2, 100 + ReadCache.DEFAULT_MAX_SOFT_CLIP));
        assertEquals(100 + ReadCache.DEFAULT_MAX_SOFT_CLIP, readCache.currentReadMinPosition());
        assertEquals(100, readCache.minCachedReadStart());
        assertEquals(2, readCache.cachedReadCount());
        assertEquals(2, readCache.cachedFragCoordGroups());
        poppedReads = collaseFragmentCoordReads(readCache.popReads());
        assertTrue(poppedReads.isEmpty());
        assertEquals(HashMultiset.create(List.of("READ_001")), readCache.auxCacheReadNames());
        assertEquals(100, readCache.minCachedReadStart());
        assertEquals(2, readCache.cachedReadCount());
        assertEquals(2, readCache.cachedFragCoordGroups());

        poppedReads = collaseFragmentCoordReads(readCache.evictAll());
        assertEquals(100 + ReadCache.DEFAULT_MAX_SOFT_CLIP, readCache.currentReadMinPosition());
        assertEquals(-1, readCache.minCachedReadStart());
        assertEquals(0, readCache.cachedReadCount());
        assertEquals(0, readCache.cachedFragCoordGroups());
        assertEquals(HashMultiset.create(List.of("READ_001", "READ_002")), poppedReads);
        assertTrue(readCache.auxCacheReadNames().isEmpty());

        poppedReads.clear();
        readCache.processRead(createRead.apply(3, 1_000));
        assertEquals(1_000, readCache.currentReadMinPosition());
        poppedReads.addAll(collaseFragmentCoordReads(readCache.popReads()));
        readCache.processRead(createRead.apply(4, 1_000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE));
        assertEquals(1_000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE, readCache.currentReadMinPosition());
        poppedReads.addAll(collaseFragmentCoordReads(readCache.popReads()));
        readCache.processRead(createRead.apply(5, 1_000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE + ReadCache.DEFAULT_MAX_SOFT_CLIP));
        assertEquals(1_000 + SINGLE_END_JITTER_COLLAPSE_DISTANCE + ReadCache.DEFAULT_MAX_SOFT_CLIP, readCache.currentReadMinPosition());
        poppedReads.addAll(collaseFragmentCoordReads(readCache.popReads()));
        assertEquals(1_000, readCache.minCachedReadStart());
        assertEquals(3, readCache.cachedReadCount());
        assertEquals(3, readCache.cachedFragCoordGroups());
        assertTrue(poppedReads.isEmpty());
        assertEquals(HashMultiset.create(List.of("READ_003", "READ_004")), readCache.auxCacheReadNames());
    }

    @Test
    public void testSingleMateUnmappedFragment()
    {
        JitterReadCache readCache = new JitterReadCache(
                new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, false, ILLUMINA));

        SAMRecord read = createSamRecord(
                "READ_001", CHR_1, 100, TEST_READ_BASES, TEST_CIGAR, CHR_1, 100, false, false, null, false, NO_CIGAR);
        read.setMateUnmappedFlag(true);

        SAMRecord mate = createSamRecord(
                "READ_001", CHR_1, 100, TEST_READ_BASES, NO_CIGAR, CHR_1, 100, false, false, null, false, TEST_CIGAR);
        mate.setReadUnmappedFlag(true);
        flipFirstInPair(mate);

        SAMRecord read2 = createSamRecord(
                "READ_002", CHR_1, 10_000, TEST_READ_BASES, TEST_CIGAR, CHR_2, 100, false, false, null, true, TEST_CIGAR);

        readCache.processRead(read);
        Multiset<String> poppedReads = collaseFragmentCoordReads(readCache.popReads());
        assertTrue(poppedReads.isEmpty());
        assertTrue(readCache.auxCacheReadNames().isEmpty());

        readCache.processRead(mate);
        poppedReads = collaseFragmentCoordReads(readCache.popReads());
        assertTrue(poppedReads.isEmpty());
        assertTrue(readCache.auxCacheReadNames().isEmpty());

        readCache.processRead(read2);
        poppedReads = collaseFragmentCoordReads(readCache.popReads());
        assertEquals(HashMultiset.create(List.of("READ_001", "READ_001")), poppedReads);
        assertTrue(readCache.auxCacheReadNames().isEmpty());
    }

    @Test
    public void testReadCacheBoundary()
    {
        ReadCache innerReadCache = new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, false, ILLUMINA);
        JitterReadCache readCache = new JitterReadCache(innerReadCache);

        SAMRecord read1 = createSamRecord(
                "READ_001", CHR_1, 1_000, TEST_READ_BASES, TEST_CIGAR, CHR_2, 100, false, false, null, true, TEST_CIGAR);

        readCache.processRead(read1);
        int expectedReadCacheBoundary = 1_000 - ReadCache.DEFAULT_MAX_SOFT_CLIP + 1;

        assertEquals(expectedReadCacheBoundary, readCache.readCacheBoundary());

        innerReadCache.evictAll();
        assertEquals(expectedReadCacheBoundary, readCache.readCacheBoundary());
    }

    @Test
    public void testDoNotGroupCoordsThatDoNotMatchOnAtLeastOneEnd()
    {
        JitterReadCache readCache = new JitterReadCache(
                new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, false, ILLUMINA));

        SAMRecord read1 = createSamRecord(
                "READ_001", CHR_1, 1_000, TEST_READ_BASES, TEST_CIGAR, CHR_2, 1_000, false, false, null, false, TEST_CIGAR);
        SAMRecord read2 = createSamRecord(
                "READ_002", CHR_1, 1_001, TEST_READ_BASES, TEST_CIGAR, CHR_2, 1_001, false, false, null, false, TEST_CIGAR);
        SAMRecord read3 = createSamRecord(
                "READ_003", CHR_1, 1_000 + ReadCache.DEFAULT_MAX_SOFT_CLIP + SINGLE_END_JITTER_COLLAPSE_DISTANCE, TEST_READ_BASES, TEST_CIGAR, CHR_2, 1_001, false, false, null, false, TEST_CIGAR);

        readCache.processRead(read1);
        readCache.processRead(read2);
        readCache.processRead(read3);

        Multiset<String> poppedReads = collaseFragmentCoordReads(readCache.popReads());

        assertEquals(HashMultiset.create(List.of("READ_001")), poppedReads);
        assertEquals(HashMultiset.create(List.of("READ_002")), readCache.auxCacheReadNames());
    }

    private static Multiset<String> collaseFragmentCoordReads(final FragmentCoordReads fragmentCoordReads)
    {
        if(fragmentCoordReads == null)
            return HashMultiset.create();

        return Stream.concat(fragmentCoordReads.SingleReads.stream().map(ReadInfo::read), fragmentCoordReads.DuplicateGroups.stream()
                        .flatMap(x -> x.reads().stream()))
                .map(SAMRecord::getReadName)
                .collect(Collectors.toCollection(HashMultiset::create));
    }
}
