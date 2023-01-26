package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.markdups.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.markdups.TestUtils.createSamRecord;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class MarkDuplicatesTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MockRefGenome mRefGenome;
    private final MarkDupsConfig mConfig;
    private final RecordWriter mWriter;
    private final PartitionDataStore mPartitionDataStore;

    public MarkDuplicatesTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mRefGenome = new MockRefGenome();
        mConfig = new MarkDupsConfig(1000, 1000, mRefGenome, true);
        mWriter = new RecordWriter(mConfig);
        mWriter.setCacheReads();
        mPartitionDataStore = new PartitionDataStore(mConfig);
    }

    @Test
    public void testSinglePartition()
    {
        ChrBaseRegion chrBaseRegion = new ChrBaseRegion(CHR_1, 1, 100000);
        ChromosomeReader chromosomeReader = new ChromosomeReader(chrBaseRegion, mConfig, mWriter, mPartitionDataStore);

        // things to test:
        // - paired and unpaired reads
        // - supplementaries and mates coming first
        // -


    }

    @Test
    public void testUmiGroup()
    {
        ChrBaseRegion chrBaseRegion = new ChrBaseRegion(CHR_1, 1, 10000);
        ChromosomeReader chromosomeReader = new ChromosomeReader(chrBaseRegion, mConfig, mWriter, mPartitionDataStore);

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

        chromosomeReader.processRead(read1);
        chromosomeReader.processRead(read2);
        chromosomeReader.flushReadPositions();

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

        chromosomeReader.processRead(mate1);
        assertEquals(3, mWriter.readsWritten().size());

        chromosomeReader.processRead(mate2);
        assertEquals(6, mWriter.readsWritten().size());

        SAMRecord supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        SAMRecord supp2 = createSamRecord(
                read2.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        chromosomeReader.processRead(supp1);
        assertEquals(6, mWriter.readsWritten().size());

        chromosomeReader.processRead(supp2);
        chromosomeReader.onChromosomeComplete();
        assertEquals(9, mWriter.readsWritten().size());
    }

    private String nextReadId(final String umiId)
    {
        return format("%s:%s", mReadIdGen.nextId(), umiId);

    }

}
