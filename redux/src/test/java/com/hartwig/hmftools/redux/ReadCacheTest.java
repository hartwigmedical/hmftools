package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.redux.duplicate.FragmentCoordReads;
import com.hartwig.hmftools.redux.duplicate.FragmentCoords;
import com.hartwig.hmftools.redux.duplicate.ReadCache;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadCacheTest
{
    @Test
    public void testBasics()
    {
        ReadCache readCache = new ReadCache(100, 100, false);

        SAMRecord read1 = createRecord(CHR_1, 50);

        readCache.processRead(read1);

        assertEquals(1, readCache.cachedReadPositionGroups());
        assertEquals(1, readCache.cachedFragCoordGroups());
        assertEquals(1, readCache.cachedReadCount());

        FragmentCoordReads fragCoordReads = readCache.popReads();
        assertNull(fragCoordReads);

        SAMRecord read2a = createRecord(CHR_1, 80);
        SAMRecord read2b = createRecord(CHR_1, 80);

        readCache.processRead(read2a);
        readCache.processRead(read2b);

        SAMRecord read3 = createRecord(CHR_1, 120);

        readCache.processRead(read3);

        assertEquals(2, readCache.cachedReadPositionGroups());
        assertEquals(4, readCache.cachedReadCount());
        assertEquals(3, readCache.cachedFragCoordGroups());

        fragCoordReads = readCache.popReads();
        assertNull(fragCoordReads);

        SAMRecord read4 = createRecord(CHR_1, 199);

        readCache.processRead(read4);

        assertEquals(2, readCache.cachedReadPositionGroups());

        // the next read triggers the first group to be processed
        SAMRecord read5 = createRecord(CHR_1, 201);

        readCache.processRead(read5);

        fragCoordReads = readCache.popReads();
        assertNotNull(fragCoordReads);

        assertEquals(1, fragCoordReads.SingleReads.size());
        assertEquals(1, fragCoordReads.DuplicateGroups.size());

        assertEquals(2, readCache.cachedReadPositionGroups());
        assertEquals(3, readCache.cachedReadCount());
        assertEquals(3, readCache.cachedFragCoordGroups());

        // next read is on a new chromosome so clears all previous
        SAMRecord read6 = createRecord(CHR_2, 50);

        readCache.processRead(read6);

        fragCoordReads = readCache.popReads();
        assertNotNull(fragCoordReads);

        assertEquals(3, fragCoordReads.SingleReads.size());
        assertEquals(0, fragCoordReads.DuplicateGroups.size());

        assertEquals(1, readCache.cachedReadPositionGroups());
        assertEquals(1, readCache.cachedReadCount());
        assertEquals(1, readCache.cachedFragCoordGroups());

        // some more duplicates in the next group
        SAMRecord read7a = createRecord(CHR_1, 150);
        SAMRecord read7b = createRecord(CHR_1, 150);
        SAMRecord read7c = createRecord(CHR_1, 150);

        readCache.processRead(read7a);
        readCache.processRead(read7b);
        readCache.processRead(read7c);

        assertEquals(2, readCache.cachedReadPositionGroups());
        assertEquals(4, readCache.cachedReadCount());
        assertEquals(2, readCache.cachedFragCoordGroups());

        fragCoordReads = readCache.evictAll();
        assertNotNull(fragCoordReads);

        assertEquals(1, fragCoordReads.SingleReads.size());
        assertEquals(1, fragCoordReads.DuplicateGroups.size());
        assertEquals(3, fragCoordReads.DuplicateGroups.get(0).totalReadCount());

        assertEquals(0, readCache.cachedReadPositionGroups());
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

        assertEquals(1, readCache.cachedReadPositionGroups());
        assertEquals(2, readCache.cachedReadCount());
        assertEquals(2, readCache.cachedFragCoordGroups());
    }

    @Test
    public void testFragmentCoordsAreOnlyPoppedOnce()
    {
        String readCigar = "12S80M51S";
        String mateCigar = "143M";

        SAMRecord read1 = createSamRecord("READ_001", CHR_1, 1, "A".repeat(143), readCigar, CHR_1, 1, false, false, null, true, mateCigar);
        SAMRecord mate1 = createSamRecord("READ_001", CHR_1, 1, "A".repeat(143), mateCigar, CHR_1, 1, true, false, null, false, readCigar);

        SAMRecord read2 = createSamRecord("READ_002", CHR_1, 1, "A".repeat(143), readCigar, CHR_1, 1, false, false, null, true, mateCigar);
        SAMRecord mate2 = createSamRecord("READ_002", CHR_1, 1, "A".repeat(143), mateCigar, CHR_1, 1, true, false, null, false, readCigar);

        FragmentCoords coords1 = FragmentCoords.fromRead(read1, false);
        FragmentCoords mateCoords1 = FragmentCoords.fromRead(mate1, false);
        FragmentCoords coords2 = FragmentCoords.fromRead(read2, false);
        FragmentCoords mateCoords2 = FragmentCoords.fromRead(mate2, false);

        assertEquals(coords1, coords2);
        assertEquals(mateCoords1, mateCoords2);
        assertNotEquals(coords1, mateCoords1);

        ReadCache readCache = new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, false);
        readCache.processRead(read1);
        readCache.processRead(mate1);
        readCache.processRead(read2);
        readCache.processRead(mate2);

        FragmentCoordReads fragmentCoordReads = readCache.evictAll();

        assertEquals(2, fragmentCoordReads.DuplicateGroups.size());
        assertEquals(0, fragmentCoordReads.SingleReads.size());

        HashMultiset<HashMultiset<String>> expectedDuplicatedGroups = HashMultiset.create(List.of(
                HashMultiset.create(List.of(read1.getSAMString(), read2.getSAMString())),
                HashMultiset.create(List.of(mate1.getSAMString(), mate2.getSAMString()))));

        HashMultiset<HashMultiset<String>> actualDuplicatedGroups = fragmentCoordReads.DuplicateGroups.stream()
                .map(g -> g.reads().stream().map(SAMRecord::getSAMString).collect(Collectors.toCollection(HashMultiset::create)))
                .collect(Collectors.toCollection(HashMultiset::create));

        assertEquals(expectedDuplicatedGroups, actualDuplicatedGroups);
    }

    @Test
    public void testReadCacheReadPositionWithinGroupBounds()
    {
        final int readLength = 143;
        final int alignmentStart = 84026397;
        final String cigar = "133M10S";
        final String readBases = "A".repeat(readLength);
        final String mateCigar = readLength + "M";

        SAMRecord read1 = createSamRecord(
                "READ_001", CHR_1, alignmentStart, readBases, cigar, CHR_2, 1, false, false, null, true, mateCigar);
        SAMRecord read2 = read1.deepCopy();
        read2.setReadName("READ_002");
        List<SAMRecord> reads = List.of(read1, read2);

        ReadCache readCache = new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, false);
        reads.forEach(readCache::processRead);
        FragmentCoordReads fragmentCoordReads = readCache.evictAll();

        assertEquals(0, fragmentCoordReads.SingleReads.size());
        assertEquals(1, fragmentCoordReads.DuplicateGroups.size());

        List<SAMRecord> duplicateGroupReads = fragmentCoordReads.DuplicateGroups.get(0).reads();
        assertEquals(2, duplicateGroupReads.size());
        assertEquals(
                Sets.newHashSet("READ_001", "READ_002"),
                duplicateGroupReads.stream().map(SAMRecord::getReadName).collect(Collectors.toCollection(Sets::newHashSet)));
    }

    private static SAMRecord createRecord(final String chromosome, final int readStart)
    {
        return createSamRecord(
                READ_ID_GEN.nextId(), chromosome, readStart, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 200, false, false, null, true, TEST_READ_CIGAR);
    }
}
