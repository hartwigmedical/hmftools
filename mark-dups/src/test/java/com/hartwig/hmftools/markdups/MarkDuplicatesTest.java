package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.markdups.TestUtils.setSecondInPair;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.markdups.common.PartitionData;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class MarkDuplicatesTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MockRefGenome mRefGenome;
    private final FileWriterCache mFileWriterCache;
    private final BamWriter mWriter;

    private final PartitionReader mPartitionReader;

    public MarkDuplicatesTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mRefGenome = new MockRefGenome();

        MarkDupsConfig config = new MarkDupsConfig(1000, 1000, mRefGenome, false, false, false);

        mFileWriterCache = new FileWriterCache(config);
        mWriter = mFileWriterCache.getBamWriter("1");
        mWriter.setCacheReads();

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
        assertEquals(1, mWriter.recordWriteCount());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate1);

        mPartitionReader.processRead(mate1);
        assertEquals(1, partitionData.resolvedFragmentStateMap().size());
        assertEquals(2, mWriter.recordWriteCount());

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mPartitionReader.postProcessRegion();
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        mPartitionReader.processRead(supp1);
        mPartitionReader.postProcessRegion();

        assertEquals(3, mWriter.recordWriteCount());
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

        assertEquals(3, mWriter.recordWriteCount());
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
        assertEquals(0, mWriter.recordWriteCount());

        // now the primaries, as candidates
        mPartitionReader.processRead(read1);
        mPartitionReader.processRead(read2);

        // mPartitionReader.flushReadPositions();

        mPartitionReader.postProcessRegion();

        assertEquals(1, partitionData.candidateDuplicatesMap().size());
        assertEquals(2, partitionData.incompleteFragmentMap().size());
        assertEquals(0, mWriter.recordWriteCount());

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
        assertEquals(0, mWriter.recordWriteCount());

        mPartitionReader.processRead(mate2);
        mPartitionReader.postProcessRegion();
        assertEquals(6, mWriter.recordWriteCount());

        assertTrue(partitionData.candidateDuplicatesMap().isEmpty());
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());
    }

    @Test
    public void testCandidateDuplicatesConsensus()
    {
        MarkDupsConfig consensusConfig = new MarkDupsConfig(
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

        assertEquals(2, mWriter.recordWriteCount());
        assertEquals(1, mWriter.recordWriteCountConsensus());

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
        assertEquals(2, mWriter.recordWriteCount());

        chrReaderConsensus.processRead(mate2);
        chrReaderConsensus.postProcessRegion();

        assertEquals(4, mWriter.recordWriteCount());
        assertEquals(2, mWriter.recordWriteCountConsensus());

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
        assertEquals(6, mWriter.recordWriteCount());
        assertEquals(3, mWriter.recordWriteCountConsensus());
    }

    @Test
    public void testExcludedRegion()
    {
        /*
        MarkDupsConfig consensusConfig = new MarkDupsConfig(
                1000, 1000, mRefGenome, false, false, true);

        PartitionDataStore partitionDataStore = new PartitionDataStore(consensusConfig);

        ChrBaseRegion chr1Region = new ChrBaseRegion(CHR_1, 1, 10000);
        ChromosomeReader chr1Reader = new ChromosomeReader(chr1Region, consensusConfig, mFileWriterCache, partitionDataStore);

        ChrBaseRegion excludedRegion = new ChrBaseRegion(CHR_2, 500, 1000);
        ChrBaseRegion chr2Region = new ChrBaseRegion(CHR_2, 1, 10000);
        ChromosomeReader chr2Reader = new ChromosomeReader(chr2Region, consensusConfig, mFileWriterCache, partitionDataStore);
        chr2Reader.setExcludedRegion(excludedRegion);

        ChrBaseRegion chr3Region = new ChrBaseRegion(CHR_3, 1, 10000);
        ChromosomeReader chr3Reader = new ChromosomeReader(chr3Region, consensusConfig, mFileWriterCache, partitionDataStore);

        int suppPos = 100;;
        int readPos = excludedRegion.start();
        int matePos = 100;

        // first fragments have supps on chr 1, primaries in the excluded regions and mates on chr 3
        String readId1 = mReadIdGen.nextId();
        String readId2 = mReadIdGen.nextId();

        SAMRecord supp1 = createSamRecord(
                readId1, CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3, matePos, false,
                true, new SupplementaryReadData(CHR_2, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord supp2 = createSamRecord(
                readId2, CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3, matePos, false,
                true, new SupplementaryReadData(CHR_2, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        chr1Reader.processRead(supp1);
        chr1Reader.processRead(supp2);
        chr1Reader.flushPendingIncompletes();

        SAMRecord read1 = createSamRecord(
                readId1, CHR_2, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), false, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                readId2, CHR_2, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), false, TEST_READ_CIGAR);

        // now the primaries, as candidates
        chr2Reader.processRead(read1);
        chr2Reader.processRead(read2);
        chr2Reader.flushReadPositions();

        PartitionData partitionData = partitionDataStore.getOrCreatePartitionData("2_0");

        assertEquals(4, mWriter.recordWriteCount());
        assertEquals(0, mWriter.recordWriteCountConsensus()); // since in excluded region

        assertEquals(2, partitionData.resolvedFragmentStateMap().size());
        assertTrue(partitionData.duplicateGroupMap().isEmpty()); // since in excluded region
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());

        SAMRecord mate1 = createSamRecord(
                readId1, CHR_3, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, readPos, true,
                false, null);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                readId2, CHR_3, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, readPos, true,
                false, null);
        setSecondInPair(mate2);

        chr3Reader.processRead(mate1);
        chr3Reader.processRead(mate2);
        chr3Reader.flushPendingIncompletes();
        assertEquals(6, mWriter.recordWriteCount());
        assertEquals(0, mWriter.recordWriteCountConsensus());

        assertTrue(partitionData.duplicateGroupMap().isEmpty());
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());

        // repeat again but with the primaries on 1, mates in the excluded region and supps on 3
        // second fragments have primaries on chr 1, mates on chr 2 in the excluded region and supps on chr 3
        mWriter.resetRecordWriteCounts();

        readPos = 100;
        matePos = excludedRegion.start();
        suppPos = 100;;

        String readId3 = mReadIdGen.nextId();
        String readId4 = mReadIdGen.nextId();

        SAMRecord read3 = createSamRecord(
                readId3, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, matePos, false, false,
                new SupplementaryReadData(CHR_3, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), false, TEST_READ_CIGAR);

        SAMRecord read4 = createSamRecord(
                readId4, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, matePos, false, false,
                new SupplementaryReadData(CHR_3, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1), false, TEST_READ_CIGAR);

        // now the primaries, as candidates
        chr1Reader.processRead(read3);
        chr1Reader.processRead(read4);
        chr1Reader.flushReadPositions();

        partitionData = partitionDataStore.getOrCreatePartitionData("1_0");

        assertEquals(2, mWriter.recordWriteCount());
        assertEquals(1, mWriter.recordWriteCountConsensus());

        assertEquals(2, partitionData.duplicateGroupMap().size());

        SAMRecord mate3 = createSamRecord(
                readId3, CHR_2, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate3);

        SAMRecord mate4 = createSamRecord(
                readId4, CHR_2, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        setSecondInPair(mate4);

        chr2Reader.processRead(mate3);
        chr2Reader.processRead(mate4);
        chr2Reader.flushPendingIncompletes();

        SAMRecord supp3 = createSamRecord(
                readId3, CHR_3, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord supp4 = createSamRecord(
                readId4, CHR_3, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        chr3Reader.processRead(supp3);
        chr3Reader.processRead(supp4);
        chr3Reader.flushPendingIncompletes();

        assertEquals(6, mWriter.recordWriteCount());
        assertEquals(3, mWriter.recordWriteCountConsensus());

        assertTrue(partitionData.duplicateGroupMap().isEmpty());
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());

        // now the primaries in the excluded region and the mates later

        String readId5 = mReadIdGen.nextId();
        String readId6 = mReadIdGen.nextId();

        readPos = 600;
        matePos = 200;

        SAMRecord read5 = createSamRecord(
                readId5, CHR_2, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3, matePos, false, false,
                null, false, TEST_READ_CIGAR);

        SAMRecord read6 = createSamRecord(
                readId6, CHR_2, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_3, matePos, false, false,
                null, false, TEST_READ_CIGAR);

        // now the primaries, as candidates
        chr2Reader.processRead(read5);
        chr2Reader.processRead(read6);
        chr2Reader.flushReadPositions();

        partitionData = partitionDataStore.getOrCreatePartitionData("2_0");

        assertEquals(2, partitionData.resolvedFragmentStateMap().size());

        assertEquals(8, mWriter.recordWriteCount());
        assertEquals(3, mWriter.recordWriteCountConsensus());

        SAMRecord mate5 = createSamRecord(
                readId5, CHR_3, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, readPos, true,
                false, null);
        setSecondInPair(mate5);

        SAMRecord mate6 = createSamRecord(
                readId6, CHR_3, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_2, readPos, true,
                false, null);
        setSecondInPair(mate6);

        chr3Reader.processRead(mate5);
        chr3Reader.processRead(mate6);
        chr3Reader.flushPendingIncompletes();

        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());

        assertEquals(10, mWriter.recordWriteCount());
        */
    }
}
