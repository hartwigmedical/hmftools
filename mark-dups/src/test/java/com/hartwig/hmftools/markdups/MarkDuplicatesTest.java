package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_CIGAR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.markdups.common.PartitionData;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class MarkDuplicatesTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MockRefGenome mRefGenome;
    private final RecordWriter mWriter;

    private final ChromosomeReader mChromosomeReader;
    private final ChromosomeReader mChrReaderUMIs;

    public MarkDuplicatesTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mRefGenome = new MockRefGenome();

        MarkDupsConfig config = new MarkDupsConfig(1000, 1000, mRefGenome, false);
        mWriter = new RecordWriter(config);
        mWriter.setCacheReads();

        ChrBaseRegion chrBaseRegion = new ChrBaseRegion(CHR_1, 1, 100000);
        mChromosomeReader = new ChromosomeReader(chrBaseRegion, config, mWriter, new PartitionDataStore(config));

        MarkDupsConfig umiConfig = new MarkDupsConfig(1000, 1000, mRefGenome, true);
        mChrReaderUMIs = new ChromosomeReader(chrBaseRegion, umiConfig, mWriter, new PartitionDataStore(umiConfig));
    }

    @Test
    public void testNonDuplicateFragment()
    {
        // test non-duplicate fragments processing reads in various orders
        int readPos = 100;
        int matePos = 200;
        int suppPos = 2000;

        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setMateNegativeStrandFlag(true);
        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mChromosomeReader.processRead(read1);
        mChromosomeReader.flushReadPositions();

        PartitionData partitionData = mChromosomeReader.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(1, partitionData.resolvedFragmentStateMap().size());
        assertEquals(1, mWriter.readsWritten().size());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate1.setFirstOfPairFlag(false);
        mate1.setSecondOfPairFlag(true);

        mChromosomeReader.processRead(mate1);
        assertEquals(1, partitionData.resolvedFragmentStateMap().size());
        assertEquals(2, mWriter.readsWritten().size());

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChromosomeReader.processRead(supp1);
        mChromosomeReader.onChromosomeComplete();
        assertEquals(3, mWriter.readsWritten().size());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());
    }

    @Test
    public void testCandidateDuplicates()
    {
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

        mChromosomeReader.processRead(supp1);
        mChromosomeReader.processRead(supp2);
        mChromosomeReader.flushPendingIncompletes();

        PartitionData partitionData = mChromosomeReader.partitionDataStore().getOrCreatePartitionData("1_1");

        assertEquals(2, partitionData.incompleteFragmentMap().size());
        assertEquals(0, mWriter.readsWritten().size());

        // now the primaries, as candidates
        mChromosomeReader.processRead(read1);
        mChromosomeReader.processRead(read2);
        mChromosomeReader.flushReadPositions();

        assertEquals(1, partitionData.candidateDuplicatesMap().size());
        assertEquals(2, partitionData.incompleteFragmentMap().size());
        assertEquals(0, mWriter.readsWritten().size());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate1.setFirstOfPairFlag(false);
        mate1.setSecondOfPairFlag(true);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate2.setFirstOfPairFlag(false);
        mate2.setSecondOfPairFlag(true);

        mChromosomeReader.processRead(mate1);
        assertEquals(0, mWriter.readsWritten().size());

        mChromosomeReader.processRead(mate2);
        mChromosomeReader.flushPendingIncompletes();
        assertEquals(6, mWriter.readsWritten().size());

        assertTrue(partitionData.candidateDuplicatesMap().isEmpty());
        assertTrue(partitionData.incompleteFragmentMap().isEmpty());
        assertTrue(partitionData.resolvedFragmentStateMap().isEmpty());
    }

    @Test
    public void testUmiGroup()
    {
        // 2 primaries in a UMI group, followed by their mates in the same partition and then their supps in a different partition
        String umidId = "TATTAT";
        int readPos = 100;
        int matePos = 200;
        int suppPos = 2000;

        SAMRecord read1 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setMateNegativeStrandFlag(true);
        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                nextReadId(umidId), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read2.setMateNegativeStrandFlag(true);
        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mChrReaderUMIs.processRead(read1);
        mChrReaderUMIs.processRead(read2);
        mChrReaderUMIs.flushReadPositions();

        PartitionData partitionData = mChrReaderUMIs.partitionDataStore().getOrCreatePartitionData("1_0");

        assertEquals(2, partitionData.umiGroupMap().size());
        assertEquals(3, mWriter.readsWritten().size());

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate1.setFirstOfPairFlag(false);
        mate1.setSecondOfPairFlag(true);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null);
        mate2.setFirstOfPairFlag(false);
        mate2.setSecondOfPairFlag(true);

        mChrReaderUMIs.processRead(mate1);
        assertEquals(3, mWriter.readsWritten().size());

        mChrReaderUMIs.processRead(mate2);
        assertEquals(6, mWriter.readsWritten().size());

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord supp2 = createSamRecord(
                read2.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        mChrReaderUMIs.processRead(supp1);
        assertEquals(6, mWriter.readsWritten().size());

        mChrReaderUMIs.processRead(supp2);
        mChrReaderUMIs.onChromosomeComplete();
        assertEquals(9, mWriter.readsWritten().size());
        assertTrue(partitionData.umiGroupMap().isEmpty());
    }

    private String nextReadId(final String umiId)
    {
        return format("%s:%s", mReadIdGen.nextId(), umiId);

    }

}
