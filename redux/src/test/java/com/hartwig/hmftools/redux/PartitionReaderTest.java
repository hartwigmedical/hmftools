package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.alignmentsToSamTag;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES_REPEAT_40;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.setSecondInPair;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class PartitionReaderTest
{
    private final ReadIdGenerator mReadIdGen;
    private final MockRefGenome mRefGenome;
    private final TestBamWriter mWriter;

    private final PartitionReader mPartitionReader;

    public PartitionReaderTest()
    {
        mReadIdGen = new ReadIdGenerator();
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES_REPEAT_40);
        mRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES_REPEAT_40.length());

        ReduxConfig umiConfig = new ReduxConfig(mRefGenome, true, false, false);

        mWriter = new TestBamWriter(umiConfig);

        mPartitionReader = new PartitionReader(umiConfig, null, mWriter, mWriter);
    }

    @Test
    public void testMultipleSupplementaries()
    {
        // one primary has multiple supplementaries and the other has a different supp
        int readPos = 100;
        int matePos = 200;
        int suppPos = 1800;

        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                null, true, TEST_READ_CIGAR);

        List<SupplementaryReadData> read1Supps = Lists.newArrayList(
                new SupplementaryReadData(CHR_1, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                new SupplementaryReadData(CHR_2, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                new SupplementaryReadData(CHR_3, suppPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1));

        read1.setAttribute(SUPPLEMENTARY_ATTRIBUTE, alignmentsToSamTag(read1Supps));
        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, readPos + 2, TEST_READ_BASES, "2S98M", CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_2, suppPos + 1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 1000));

        mPartitionReader.processRead(read1);
        mPartitionReader.processRead(read2);
        mPartitionReader.flushReadPositions();

        SAMRecord mate1 = createSamRecord(
                read1.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, null, false, TEST_READ_CIGAR);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                read2.getReadName(), CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos + 2, true,
                false, null, false, TEST_READ_CIGAR);
        setSecondInPair(mate2);

        mPartitionReader.processRead(mate1);
        mPartitionReader.processRead(mate2);

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1001, 2000));

        SAMRecord read1Supp1 = createSamRecord(
                read1.getReadName(), CHR_1, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReader.processRead(read1Supp1);
        mPartitionReader.postProcessRegion();

        SAMRecord read1Supp2 = createSamRecord(
                read1.getReadName(), CHR_2, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        SAMRecord read2Supp1 = createSamRecord(
                read2.getReadName(), CHR_2, suppPos + 1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos + 2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_2, 1001, 2000));

        mPartitionReader.processRead(read1Supp2);
        mPartitionReader.processRead(read2Supp1);
        mPartitionReader.postProcessRegion();

        SAMRecord read1Supp3 = createSamRecord(
                read1.getReadName(), CHR_3, suppPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_3, 1001, 2000));

        mPartitionReader.processRead(read1Supp3);
        mPartitionReader.postProcessRegion();

        assertEquals(8, mWriter.nonConsensusWriteCount());
        assertEquals(1, mWriter.consensusWriteCount()); // since the supplementaries aren't mapped to the same location
    }

    @Test
    public void testInconsistentSupplementaries()
    {
        // supplementaries arriving out of order and mapped to different locations
        int supp2Pos1 = 950;
        int readPos = 1950;
        int supp2Pos2 = 2950;
        int matePos = 3500;
        int supp1Pos1 = 4500;
        int supp1Pos2 = 5500;

        String readId1 = mReadIdGen.nextId();
        String readId2 = mReadIdGen.nextId();

        SAMRecord read2Supp1 = createSamRecord(
                readId2, CHR_1, supp2Pos1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, 1, 999));
        mPartitionReader.processRead(read2Supp1);
        mPartitionReader.postProcessRegion();

        SAMRecord read1 = createSamRecord(
                readId1, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, supp1Pos2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        read1.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        SAMRecord read2 = createSamRecord(
                readId2, CHR_1, readPos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false, false,
                new SupplementaryReadData(CHR_1, supp2Pos1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        read2.setAttribute(MATE_CIGAR_ATTRIBUTE, TEST_READ_CIGAR);

        int regionStart = 1000;
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReader.processRead(read1);
        mPartitionReader.processRead(read2);

        mPartitionReader.postProcessRegion();

        SAMRecord read2Supp2 = createSamRecord(
                readId2, CHR_1, supp2Pos2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);
        setSecondInPair(read2Supp2);

        regionStart += 1000;
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReader.processRead(read2Supp2);
        mPartitionReader.postProcessRegion();

        SAMRecord mate1 = createSamRecord(
                readId1, CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, supp1Pos1, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                false, TEST_READ_CIGAR);
        setSecondInPair(mate1);

        SAMRecord mate2 = createSamRecord(
                readId2, CHR_1, matePos, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, true,
                false, new SupplementaryReadData(CHR_1, supp2Pos2, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                false, TEST_READ_CIGAR);
        setSecondInPair(mate2);

        regionStart += 1000;
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReader.processRead(mate1);
        mPartitionReader.processRead(mate2);
        mPartitionReader.postProcessRegion();

        SAMRecord read1Supp1 = createSamRecord(
                read1.getReadName(), CHR_1, supp1Pos2, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, readPos, false,
                true, new SupplementaryReadData(CHR_1, matePos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);
        setSecondInPair(read1Supp1);

        SAMRecord read1Supp2 = createSamRecord(
                read1.getReadName(), CHR_1, supp1Pos1, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, matePos, false,
                true, new SupplementaryReadData(CHR_1, readPos, SUPP_POS_STRAND, TEST_READ_CIGAR, 1),
                true, TEST_READ_CIGAR);

        // note last 2 regions are out of order
        regionStart = 5000;
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReader.processRead(read1Supp1);
        mPartitionReader.postProcessRegion();

        regionStart = 4000;
        mPartitionReader.setupRegion(new ChrBaseRegion(CHR_1, regionStart, regionStart + 999));
        mPartitionReader.processRead(read1Supp2);
        mPartitionReader.postProcessRegion();

        assertEquals(8, mWriter.nonConsensusWriteCount());
        assertEquals(2, mWriter.consensusWriteCount());
    }
}
