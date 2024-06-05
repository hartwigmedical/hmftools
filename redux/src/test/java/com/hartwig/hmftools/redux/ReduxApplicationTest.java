package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES_REPEAT_40;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_MAX_NON_OVERLAPPING_BASES;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_MIN_HIGH_DEPTH;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.redux.common.HighDepthRegion;
import com.hartwig.hmftools.redux.common.PartitionData;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReduxApplicationTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MockRefGenome mRefGenome;
    // private final FileWriterCache mFileWriterCache;
    private final TestBamWriter mWriter;

    private final PartitionReader mPartitionReader;

    public ReduxApplicationTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES_REPEAT_40);
        mRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES_REPEAT_40.length());

        ReduxConfig config = new ReduxConfig(1000, 1000, mRefGenome, false, false, false);

        // mFileWriterCache = new FileWriterCache(config);
        // mWriter = mFileWriterCache.getPartitionBamWriter("1");
        mWriter = new TestBamWriter(config);

        mPartitionReader = new PartitionReader(config, null, mWriter, new PartitionDataStore(config));
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

        PartitionData partitionData = mPartitionReader.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(1, partitionData.resolvedFragmentStateMap().size());
        assertEquals(1, mWriter.nonConsensusWriteCount());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate1);

        mPartitionReader.processRead(mate1);
        assertEquals(1, partitionData.resolvedFragmentStateMap().size());
        assertEquals(2, mWriter.nonConsensusWriteCount());

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mPartitionReader.postProcessRegion();
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        mPartitionReader.processRead(supp1);
        mPartitionReader.postProcessRegion();

        assertEquals(3, mWriter.nonConsensusWriteCount());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());
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
                false, null);

        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mPartitionReader.processRead(read1);

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mPartitionReader.processRead(supp1);

        mPartitionReader.postProcessRegion();

        PartitionData partitionData = mPartitionReader.partitionDataStore().getOrCreatePartitionData("1_0");
        partitionData.incompleteFragmentMap().containsKey(supp1.getReadName());

        partitionData.incompleteFragmentMap().isEmpty();

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));
        setSecondInPair(mate1);

        mPartitionReader.processRead(mate1);
        mPartitionReader.postProcessRegion();

        assertEquals(3, mWriter.nonConsensusWriteCount());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());
    }

    @Test
    public void testCandidateDuplicates()
    {
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        // no use of mate CIGAR so rely on putting reads together to determine duplicate status
        int suppPos = 500;
        int readPos = 1500;
        int matePos = 2500;

        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setMateNegativeStrandFlag(true);

        SAMRecord read2 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read2.setMateNegativeStrandFlag(true);

        // send through supplementaries first
        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord supp2 = createSamRecord(
                read2.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mPartitionReader.processRead(supp1);
        mPartitionReader.processRead(supp2);

        // mPartitionReader.flushPendingIncompletes();

        mPartitionReader.postProcessRegion();

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        PartitionData partitionData = mPartitionReader.partitionDataStore().getOrCreatePartitionData("1_1");

        assertEquals(2, partitionData.incompleteFragmentMap().size());
        assertEquals(0, mWriter.nonConsensusWriteCount());

        // now the primaries, as candidates
        mPartitionReader.processRead(read1);
        mPartitionReader.processRead(read2);

        // mPartitionReader.flushReadPositions();

        mPartitionReader.postProcessRegion();

        assertEquals(1, partitionData.candidateDuplicatesMap().size());
        assertEquals(2, partitionData.incompleteFragmentMap().size());
        assertEquals(0, mWriter.nonConsensusWriteCount());

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 2001, 3000));

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate2);

        mPartitionReader.processRead(mate1);
        assertEquals(0, mWriter.nonConsensusWriteCount());

        mPartitionReader.processRead(mate2);
        mPartitionReader.postProcessRegion();
        assertEquals(6, mWriter.nonConsensusWriteCount());

        assertTrue(partitionData.candidateDuplicatesMap().isEmpty());
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());
    }

    @Test
    public void testCandidateDuplicatesConsensus()
    {
        ReduxConfig consensusConfig = new ReduxConfig(
                1000, 1000, mRefGenome, false, false, true);

        PartitionReader chrReaderConsensus = new PartitionReader(
                consensusConfig, null, mWriter, new PartitionDataStore(consensusConfig));

        chrReaderConsensus.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        int readPos = 100;
        int matePos = 1500;
        int suppPos = 2500;

        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        // now the primaries, as candidates
        chrReaderConsensus.processRead(read1);
        chrReaderConsensus.processRead(read2);
        chrReaderConsensus.postProcessRegion();

        chrReaderConsensus.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        PartitionData partitionData = chrReaderConsensus.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(2, mWriter.nonConsensusWriteCount());
        assertEquals(1, mWriter.consensusWriteCount());

        assertEquals(2, partitionData.duplicateGroupMap().size());
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());
        // assertEquals(0, mWriter.recordWriteCount());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate2);

        chrReaderConsensus.processRead(mate1);
        assertEquals(2, mWriter.nonConsensusWriteCount());

        chrReaderConsensus.processRead(mate2);
        chrReaderConsensus.postProcessRegion();

        assertEquals(4, mWriter.nonConsensusWriteCount());
        assertEquals(2, mWriter.consensusWriteCount());

        // send through supplementaries first
        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord supp2 = createSamRecord(
                read2.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        chrReaderConsensus.setupRegion(new ChrBaseRegion(CHR_1, 2001, 3000));

        chrReaderConsensus.processRead(supp1);
        chrReaderConsensus.processRead(supp2);
        chrReaderConsensus.postProcessRegion();

        assertTrue(partitionData.duplicateGroupMap().isEmpty());
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());
        assertEquals(6, mWriter.nonConsensusWriteCount());
        assertEquals(3, mWriter.consensusWriteCount());
    }

    @Test
    public void testUnmapRegionReads()
    {
        ReduxConfig config = new ReduxConfig(
                2000, 1000, mRefGenome, false, false, false);

        HighDepthRegion region1 = new HighDepthRegion(550, 650, UNMAP_MIN_HIGH_DEPTH + 1);
        HighDepthRegion region2 = new HighDepthRegion(900, 1300, UNMAP_MIN_HIGH_DEPTH + 1);
        config.UnmapRegions.addRegion(CHR_1, region1);
        config.UnmapRegions.addRegion(CHR_1, region2);

        PartitionReader partitionReader = new PartitionReader(config, null, mWriter, new PartitionDataStore(config));

        partitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1999));

        // first fragments have supps on chr 1, primaries in the excluded regions and mates on chr 3

        // a read not overlapping enough with a region to be unmapped
        // but supp is unmapped
        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1,
                550 - UNMAP_MAX_NON_OVERLAPPING_BASES - 2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3,
                800, false, false,
                new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                false, TEST_READ_CIGAR);

        partitionReader.processRead(read1);

        assertEquals(null, read1.getAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertFalse(read1.getReadUnmappedFlag());

        // Read overlaps so that the non overlapping bases of the aligned part does not exceed `ReadUnmapper.MAX_NON_OVERLAPPING_BASES`.
        SAMRecord read2 = createSamRecord(
                mReadIdGen.nextId(), CHR_1,
                551 + UNMAP_MAX_NON_OVERLAPPING_BASES - 1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3,
                800, false, false,
                null, false, TEST_READ_CIGAR);

        partitionReader.processRead(read2);

        assertTrue(read2.getReadUnmappedFlag());

        // mate is unmapped
        SAMRecord read3 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1,
                550 - UNMAP_MAX_NON_OVERLAPPING_BASES + 1, false, false,
                null, false, TEST_READ_CIGAR);

        partitionReader.processRead(read3);

        assertFalse(read3.getReadUnmappedFlag());
        assertTrue(read3.getMateUnmappedFlag());

        // unmapped read then has mate unmapped
        SAMRecord read4 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 1000, TEST_READ_BASES, NO_CIGAR, CHR_1,
                1000, false, false,
                null, false, TEST_READ_CIGAR);

        read4.setReadUnmappedFlag(true);

        partitionReader.processRead(read4);

        assertTrue(read4.getReadUnmappedFlag());
        assertTrue(read4.getMateUnmappedFlag());
    }
}
