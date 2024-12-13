package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.PartitionThread.splitRegionsIntoPartitions;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES_REPEAT_40;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createPartitionRead;
import static com.hartwig.hmftools.redux.TestUtils.createTestConfig;
import static com.hartwig.hmftools.redux.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_MAX_NON_OVERLAPPING_BASES;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_MIN_HIGH_DEPTH;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.UnmappingRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.unmap.ReadUnmapper;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReduxApplicationTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MockRefGenome mRefGenome;
    private final TestBamWriter mWriter;

    private final PartitionReader mPartitionReader;

    public ReduxApplicationTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES_REPEAT_40);
        mRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES_REPEAT_40.length());

        ReduxConfig config = createTestConfig();

        mWriter = new TestBamWriter(config);

        mPartitionReader = createPartitionRead(config, mWriter);
    }

    @Test
    public void testNonDuplicateFragment()
    {
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        // test non-duplicate fragments processing reads in various orders
        int readPos = 100;
        int matePos = 200;
        int suppPos = 2000;

        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), true, TEST_READ_CIGAR);

        mPartitionReader.processRead(read1);
        mPartitionReader.flushReadPositions();

        assertEquals(1, mWriter.nonConsensusWriteCount());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null, false, TEST_READ_CIGAR);
        setSecondInPair(mate1);

        mPartitionReader.processRead(mate1);

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReader.postProcessRegion();
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        mPartitionReader.processRead(supp1);
        mPartitionReader.postProcessRegion();

        assertEquals(3, mWriter.nonConsensusWriteCount());
    }

    @Test
    public void testSupplementaryFragments()
    {
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        // first a supplementary linked to a distant mate but in the same partition as the primary
        int readPos = 100;
        int matePos = 1800;
        int suppPos = 200;

        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                false, null, true, TEST_READ_CIGAR);

        mPartitionReader.processRead(read1);

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReader.processRead(supp1);

        mPartitionReader.postProcessRegion();

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                false, TEST_READ_CIGAR);
        setSecondInPair(mate1);

        mPartitionReader.processRead(mate1);
        mPartitionReader.postProcessRegion();

        assertEquals(3, mWriter.nonConsensusWriteCount());
    }

    @Test
    public void testUnmapRegionReads()
    {
        Map<String, List<UnmappingRegion>> chrLocationsMap = Maps.newHashMap();
        UnmappingRegion region1 = new UnmappingRegion(550, 650, UNMAP_MIN_HIGH_DEPTH + 1);
        UnmappingRegion region2 = new UnmappingRegion(900, 1300, UNMAP_MIN_HIGH_DEPTH + 1);
        chrLocationsMap.put(CHR_1, Lists.newArrayList(region1, region2));

        ReadUnmapper readUnmapper = new ReadUnmapper(chrLocationsMap);
        ReduxConfig config = new ReduxConfig(mRefGenome, false, false, false, readUnmapper);

        PartitionReader partitionReader = createPartitionRead(config, mWriter);

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1999));

        // first fragments have supps on chr 1, primaries in the excluded regions and mates on chr 3

        // a read not overlapping enough with a region to be unmapped but supp is unmapped
        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1,
                550 - UNMAP_MAX_NON_OVERLAPPING_BASES - 2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3,
                800, false, false,
                new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                false, TEST_READ_CIGAR);

        partitionReader.processRead(read1);

        assertEquals(null, read1.getAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertFalse(read1.getReadUnmappedFlag());

        // Read overlaps so that the non overlapping bases of the aligned part does not exceed `ReadUnmapper.MAX_NON_OVERLAPPING_BASES`.
        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1,
                551 + UNMAP_MAX_NON_OVERLAPPING_BASES - 1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3,
                800, false, false,
                null, false, TEST_READ_CIGAR);

        partitionReader.processRead(read2);

        assertTrue(read2.getReadUnmappedFlag());

        // mate is unmapped
        SAMRecord read3 = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                550 - UNMAP_MAX_NON_OVERLAPPING_BASES + 1, false, false,
                null, false, TEST_READ_CIGAR);

        partitionReader.processRead(read3);

        assertFalse(read3.getReadUnmappedFlag());
        assertTrue(read3.getMateUnmappedFlag());

        // unmapped read then has mate unmapped
        SAMRecord read4 = SamRecordTestUtils.createSamRecord(
                mReadIdGen.nextId(), CHR_1, 1000, TEST_READ_BASES, NO_CIGAR, CHR_1,
                1000, false, false,
                null, false, TEST_READ_CIGAR);

        read4.setReadUnmappedFlag(true);

        partitionReader.processRead(read4);

        assertTrue(read4.getReadUnmappedFlag());
        assertTrue(read4.getMateUnmappedFlag());
    }

    @Test
    public void testRegionPartitionSplits()
    {
        SpecificRegions specificRegions = new SpecificRegions();

        specificRegions.addRegion(new ChrBaseRegion(CHR_1, 1, 1000000));

        List<List<ChrBaseRegion>> partitionRegions = splitRegionsIntoPartitions(specificRegions, 4, V37, null);
        assertEquals(4, partitionRegions.size());

        ChrBaseRegion region = partitionRegions.get(0).get(0);
        assertEquals(1, region.start());
        assertEquals(250000, region.end());

        region = partitionRegions.get(1).get(0);
        assertEquals(250001, region.start());
        assertEquals(500000, region.end());

        region = partitionRegions.get(2).get(0);
        assertEquals(500001, region.start());
        assertEquals(750000, region.end());

        region = partitionRegions.get(3).get(0);
        assertEquals(750001, region.start());
        assertEquals(1000000, region.end());

        specificRegions = new SpecificRegions();
        specificRegions.addRegion(new ChrBaseRegion(CHR_1, 1, 10000));
        specificRegions.addRegion(new ChrBaseRegion(CHR_2, 1, 10000));
        specificRegions.addRegion(new ChrBaseRegion(CHR_3, 1, 10000));

        partitionRegions = splitRegionsIntoPartitions(specificRegions, 1, V37, null);
        assertEquals(1, partitionRegions.size());

        region = partitionRegions.get(0).get(0);
        assertEquals(specificRegions.Regions.get(0), partitionRegions.get(0).get(0));
        assertEquals(specificRegions.Regions.get(1), partitionRegions.get(0).get(1));
        assertEquals(specificRegions.Regions.get(2), partitionRegions.get(0).get(2));

        // TODO: test splits across chromosome and region boundaries, and that lengths are roughly equal
    }
}
