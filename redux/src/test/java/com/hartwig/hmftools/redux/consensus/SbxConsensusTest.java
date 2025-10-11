package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractConsensusType;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_ADJACENT_1_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_ADJACENT_2_3_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_READ_INDEX_TAG;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.extractDuplexBaseIndex;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildBaseQuals;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.cloneSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.createConsensusRead;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MATCH;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_SOFTCLIP;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.DUPLEX_NO_CONSENSUS_QUAL;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.SIMPLEX_NO_CONSENSUS_QUAL;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.ReduxConstants;

import org.junit.After;
import org.junit.Ignore;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SbxConsensusTest
{
    private final MockRefGenome mRefGenome;
    private final ConsensusReads mConsensusReads;

    //                                                 20        30        40        50
    //                                        12345678901234567890123456789012345678901234567890
    private static final String REF_BASES = "XAAAAACCCCCGGGGGTTTTTACGTACGTACGTAACCGGTTAACCGGTT";

    public SbxConsensusTest()
    {
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        mRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES.length());
        mConsensusReads = new ConsensusReads(mRefGenome, SequencingType.SBX, new ConsensusStatistics());
        mConsensusReads.setChromosomeLength(REF_BASES.length());

        ReduxConfig.SEQUENCING_TYPE = SequencingType.SBX;
        ReduxConfig.RunChecks = true;
    }

    @After
    public void resetSequencingType()
    {
        ReduxConfig.SEQUENCING_TYPE = SequencingType.ILLUMINA;
        ReduxConfig.RunChecks = false;
    }

    @Test
    public void testSbxSimplexConsensus()
    {
        List<SAMRecord> reads = Lists.newArrayList();

        int position = 1;
        String refBases = REF_BASES.substring(position, 11);

        // reads have a single disagreeing base but find > 50% in agreement
        byte[] baseQuals = buildBaseQuals(refBases.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read1 = createSamRecord(refBases, position, baseQuals);
        reads.add(read1);

        String readBases = "T" + refBases.substring(1);
        baseQuals = buildBaseQuals(readBases.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read2 = createSamRecord(readBases, position, baseQuals);
        reads.add(read2);

        baseQuals = buildBaseQuals(readBases.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read3 = createSamRecord(readBases, position, baseQuals);
        reads.add(read3);

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(ALIGNMENT_ONLY, readInfo.Outcome);
        assertEquals(readBases, readInfo.ConsensusRead.getReadString());
        assertEquals("10M", readInfo.ConsensusRead.getCigarString());

        for(int i = 0; i < readInfo.ConsensusRead.getBaseQualities().length; ++i)
        {
            assertEquals(RAW_SIMPLEX_QUAL, readInfo.ConsensusRead.getBaseQualities()[i]);
        }

        assertEquals(position, readInfo.ConsensusRead.getAlignmentStart());

        // test again with no bases agreeing
        readBases = "G" + refBases.substring(1);

        read3 = createSamRecord(readBases, position, baseQuals);
        reads.set(2, read3);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(refBases, readInfo.ConsensusRead.getReadString());
        assertEquals(SIMPLEX_NO_CONSENSUS_QUAL, readInfo.ConsensusRead.getBaseQualities()[0]);
    }

    @Test
    public void testSbxDuplexConsensus()
    {
        List<SAMRecord> reads = Lists.newArrayList();

        int position = 1;
        String refBases = REF_BASES.substring(position, 11);

        String readBases = "T" + refBases.substring(1);

        // reads all agree
        byte[] baseQuals = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);
        SAMRecord read1 = createSamRecord(readBases, position, baseQuals);
        reads.add(read1);

        baseQuals = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);
        SAMRecord read2 = createSamRecord(readBases, position, baseQuals);
        reads.add(read2);

        baseQuals = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);
        SAMRecord read3 = createSamRecord(readBases, position, baseQuals);
        reads.add(read3);

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(readBases, readInfo.ConsensusRead.getReadString());
        assertEquals("10M", readInfo.ConsensusRead.getCigarString());

        for(int i = 0; i < readInfo.ConsensusRead.getBaseQualities().length; ++i)
        {
            assertEquals(RAW_DUPLEX_QUAL, readInfo.ConsensusRead.getBaseQualities()[i]);
        }

        assertEquals(position, readInfo.ConsensusRead.getAlignmentStart());

        // test again with 2 of 3 agreeing, so a consensus is used
        String readBases3 = "G" + refBases.substring(1);

        read3 = createSamRecord(readBases3, position, baseQuals);
        reads.set(2, read3);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(readBases, readInfo.ConsensusRead.getReadString());
        assertEquals(RAW_DUPLEX_QUAL, readInfo.ConsensusRead.getBaseQualities()[0]);

        // test that low qual reads are considered for the required percentage but not consensus itself
        baseQuals = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);
        baseQuals[0] = SBX_DUPLEX_MISMATCH_QUAL;
        SAMRecord read4 = createSamRecord(readBases, position, baseQuals);
        reads.add(read4);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(refBases, readInfo.ConsensusRead.getReadString());
        assertEquals(DUPLEX_NO_CONSENSUS_QUAL, readInfo.ConsensusRead.getBaseQualities()[0]);

        // a simplex read is not enough for consensus again
        baseQuals = buildBaseQuals(readBases.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read5 = createSamRecord(readBases, position, baseQuals);
        reads.add(read5);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(refBases, readInfo.ConsensusRead.getReadString());
        assertEquals(DUPLEX_NO_CONSENSUS_QUAL, readInfo.ConsensusRead.getBaseQualities()[0]);

        baseQuals = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);
        SAMRecord read6 = createSamRecord(readBases, position, baseQuals);
        reads.add(read6);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(readBases, readInfo.ConsensusRead.getReadString());
        assertEquals(RAW_DUPLEX_QUAL, readInfo.ConsensusRead.getBaseQualities()[0]);
    }

    @Test
    public void testSbxSimpleIndelConsensus()
    {
        List<SAMRecord> reads = Lists.newArrayList();

        int position = 1;

        String readBasesStart = "AACC";
        String readBasesEnd = "GGTT";

        // test 1: reads differ by an inserted T only, favours M over I since all counts are equal
        String readBases1 = readBasesStart + "TTTT" + readBasesEnd;
        String readCigar1 = "4M1I7M";

        byte[] readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_DUPLEX_QUAL);
        SAMRecord read1 = createSamRecord(readBases1, position, readBaseQuals1, readCigar1);
        reads.add(read1);

        String readBases2 = readBasesStart + "TTT" + readBasesEnd;
        String readCigar2 = "11M";

        byte[] readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_DUPLEX_QUAL);
        SAMRecord read2 = createSamRecord(readBases2, position, readBaseQuals2, readCigar2);
        reads.add(read2);

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readBases2, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigar2, readInfo.ConsensusRead.getCigarString());
        assertEquals(position, readInfo.ConsensusRead.getAlignmentStart());

        /*
        for(int i = 0; i < readInfo.ConsensusRead.getBaseQualities().length; ++i)
        {
            assertEquals(RAW_SIMPLEX_QUAL, readInfo.ConsensusRead.getBaseQualities()[i]);
        }
        */

        // test 2: favour ref over insert from 2 duplex reads, ignoring simplex read
        reads.clear();

        readBases1 = readBasesStart + "TTTT" + readBasesEnd;
        readCigar1 = "4M1I7M";
        readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_DUPLEX_QUAL);
        read1 = createSamRecord(readBases1, position, readBaseQuals1, readCigar1);
        reads.add(read1);

        readBases2 = readBases1;
        readCigar2 = readCigar1;
        readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_DUPLEX_QUAL);
        readBaseQuals2[4] = RAW_SIMPLEX_QUAL;
        read2 = createSamRecord(readBases2, position, readBaseQuals2, readCigar2);
        reads.add(read2);

        String readBases3 = readBasesStart + "TTT" + readBasesEnd;
        String readCigar3 = "11M";
        byte[] readBaseQuals3 = buildBaseQuals(readBases2.length(), RAW_DUPLEX_QUAL);
        SAMRecord read3 = createSamRecord(readBases3, position, readBaseQuals3, readCigar3);
        reads.add(read3);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readBases3, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigar3, readInfo.ConsensusRead.getCigarString());
        assertEquals(position, readInfo.ConsensusRead.getAlignmentStart());

        // test 3: both reads have the same cigar, so reverts to standard qual-checking logic
        reads.clear();

        readBases1 = readBasesStart + "TTTT" + readBasesEnd;
        readCigar1 = "4M1I7M";
        readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_DUPLEX_QUAL);
        read1 = createSamRecord(readBases1, position, readBaseQuals1, readCigar1);
        reads.add(read1);

        readBases2 = readBasesStart + "GTTT" + readBasesEnd;
        readCigar2 = readCigar1;
        readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_DUPLEX_QUAL);
        readBaseQuals2[4] = RAW_DUPLEX_MISMATCH_QUAL;
        read2 = createSamRecord(readBases2, position, readBaseQuals2, readCigar2);
        reads.add(read2);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MATCH, readInfo.Outcome);
        String readBasesNet = readBasesStart + "ATTT" + readBasesEnd; // using ref since no agreement
        assertEquals(readBasesNet, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigar2, readInfo.ConsensusRead.getCigarString());
        assertEquals(DUPLEX_NO_CONSENSUS_QUAL, readInfo.ConsensusRead.getBaseQualities()[4]);

        // test 5:
        reads.clear();

        readBases1 = readBasesStart + "TTTT" + readBasesEnd;
        readCigar1 = "4M1D8M";
        readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_DUPLEX_QUAL);
        read1 = createSamRecord(readBases1, position, readBaseQuals1, readCigar1);
        reads.add(read1);

        readBases2 = readBasesStart + "TTT" + readBasesEnd;
        readCigar2 = "4M2D7M";
        readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_DUPLEX_QUAL);
        read2 = createSamRecord(readBases2, position, readBaseQuals2, readCigar2);
        reads.add(read2);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readBases1, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigar1, readInfo.ConsensusRead.getCigarString());

        // test 6:
        reads.clear();

        readBases1 = readBasesStart + "TTTT" + readBasesEnd;
        readCigar1 = "4M1I7M";
        readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_DUPLEX_QUAL);
        read1 = createSamRecord(readBases1, position, readBaseQuals1, readCigar1);
        reads.add(read1);

        readBases2 = readBasesStart + "TTTTT" + readBasesEnd;
        readCigar2 = "4M2I7M";
        readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_DUPLEX_QUAL);
        read2 = createSamRecord(readBases2, position, readBaseQuals2, readCigar2);
        reads.add(read2);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readBases1, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigar1, readInfo.ConsensusRead.getCigarString());
    }

    @Test
    public void testSbxIndelConsensus2()
    {
        List<SAMRecord> reads = Lists.newArrayList();

        int position = 1;

        // 1:   AAC GGT AAC GGT
        //      MMM MMMDMMM MMM

        // 2:   AACCGGT AAC GGT
        //      MMMIMMMDMMM MMM

        // 3:   AACCGGTTAACCGGT
        //      MMMIMMMMMMMIMMM

        // net: AACCGGT AAC GGT
        //      MMMIMMMDMMM MMM

        String readBases1 = "AACGGTAACGGT";
        String readCigar1 = "6M1D6M";

        String readBases2 = "AACCGGTAACGGT";
        String readCigar2 = "3M1I3M1D6M";

        String readBases3 = "AACCGGTTAACCGGT";
        String readCigar3 = "3M1I7M1I3M";

        String readBasesNet = "AACCGGTAACGGT";
        String readCigarNet = "3M1I3M1D6M";

        byte[] readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read1 = createSamRecord(readBases1, position, readBaseQuals1, readCigar1);
        reads.add(read1);

        byte[] readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read2 = createSamRecord(readBases2, position, readBaseQuals2, readCigar2);
        reads.add(read2);

        byte[] readBaseQuals3 = buildBaseQuals(readBases3.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read3 = createSamRecord(readBases3, position, readBaseQuals3, readCigar3);
        reads.add(read3);

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readBasesNet, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigarNet, readInfo.ConsensusRead.getCigarString());

        //                      duplex->
        // 1:   AAC   CGGTTAA A C      CGGTTA
        // 2:   AAC   CGGTTAA   C CCCC CGGTTA
        // 3:   AAC C CGGTTAA   C       GGTTA
        //       duplex->
        // net: AAC   CGGTTAA   C       GGTTA

        //            0123456789     0     1234567
        readBases1 = "AACCGGTTAA" + "A" + "CCGGTTA";
        readCigar1 = "10M1I7M";

        //            01234567890     1234     567980
        readBases2 = "AACCGGTTAAC" + "CCCC" + "CGGTTA";
        readCigar2 = "11M4I6M";

        //            012     3     45678901   D 23456
        readBases3 = "AAC" + "C" + "CGGTTAAC" + "GGTTA";
        readCigar3 = "3M1I8M1D5M";

        //              01234567890  D  12345
        readBasesNet = "AACCGGTTAAC" + "GGTTA";
        readCigarNet = "11M1D5M";

        reads.clear();

        readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_SIMPLEX_QUAL);
        setBaseQuals(readBaseQuals1, 16, -1, RAW_DUPLEX_QUAL);
        read1 = createSamRecord(readBases1, position, readBaseQuals1, readCigar1);
        reads.add(read1);

        readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_SIMPLEX_QUAL);
        setBaseQuals(readBaseQuals2, 16, -1, RAW_DUPLEX_QUAL);
        read2 = createSamRecord(readBases2, position, readBaseQuals2, readCigar2);
        reads.add(read2);

        readBaseQuals3 = buildBaseQuals(readBases3.length(), RAW_SIMPLEX_QUAL);
        setBaseQuals(readBaseQuals3, 8, -1, RAW_DUPLEX_QUAL);
        read3 = createSamRecord(readBases3, position, readBaseQuals3, readCigar3);
        reads.add(read3);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readBasesNet, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigarNet, readInfo.ConsensusRead.getCigarString());
    }

    private static void setBaseQuals(final byte[] quals, int startIndex, int endIndex, byte qual)
    {
        if(endIndex < 0)
            endIndex = quals.length - 1;

        for(int i = startIndex; i <= endIndex; ++i)
        {
            quals[i] = qual;
        }
    }

    @Test
    public void testSbxIndelConsensusWithSoftClips()
    {
        // requires ability to have different unclipped and/or aligned start locations
        List<SAMRecord> reads = Lists.newArrayList();

        // 1:   AACA GGGTT AA
        //      SSSS MMIMM MS

        // 2:     CC GGTT TACC
        //        SM MMMM SSSS

        // net: AACC GGTT AACC
        //      SSSM MMMM MSSS

        String readBases1 = "AACAGGGTTAA";
        String readCigar1 = "4S2M1I3M1S";
        int readStart1 = 10;

        String readBases2 = "CCGGTTTACC";
        String readCigar2 = "1S5M4S";
        int readStart2 = 9;

        byte[] readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read1 = createSamRecord(readBases1, readStart1, readBaseQuals1, readCigar1);
        reads.add(read1);

        byte[] readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read2 = createSamRecord(readBases2, readStart2, readBaseQuals2, readCigar2);
        reads.add(read2);

        String readBasesNet = "AACCGGTTAACC"; // chooses ref since both are simplex
        String readCigarNet = "3S6M3S";

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_SOFTCLIP, readInfo.Outcome);
        assertEquals(readStart2, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals(readBases2, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigar2, readInfo.ConsensusRead.getCigarString());

        // now with 3+ reads
        SAMRecord read3 = cloneSamRecord(read1, READ_ID_GEN.nextId());
        SAMRecord read4 = cloneSamRecord(read2, READ_ID_GEN.nextId());
        reads.add(read3);
        reads.add(read4);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readStart2, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals(readBasesNet, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigarNet, readInfo.ConsensusRead.getCigarString());

        // test 2: test in reverse direction
        //reads.clear();
        reads.forEach(x -> x.setReadNegativeStrandFlag(true));

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readStart2, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals(readBasesNet, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigarNet, readInfo.ConsensusRead.getCigarString());

        reads.clear();

        // soft-clips vs indels - indels are ignored from ref base advancement

        // 1x2  : AAACC GGTT AAACC
        //        MMIMM MMMM MMIMM

        // 3:     TTTT GGTT GGGG
        //        SSSS MMMM SSSS

        // net: to match read 1

        readBases1 = "AAACCGGTTAAACC";
        readCigar1 = "2M1I8M1I2M";
        readStart1 = 6;

        readBases2 = "TTTTGGTTGGGG";
        readCigar2 = "4S4M4S";
        readStart2 = 10;

        readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_SIMPLEX_QUAL);
        read1 = createSamRecord(readBases1, readStart1, readBaseQuals1, readCigar1);
        reads.add(read1);
        reads.add(cloneSamRecord(read1, READ_ID_GEN.nextId()));

        readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_SIMPLEX_QUAL);
        read2 = createSamRecord(readBases2, readStart2, readBaseQuals2, readCigar2);
        reads.add(read2);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readStart1, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals(readBases1, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigar1, readInfo.ConsensusRead.getCigarString());

        // and in reverse
        reads.forEach(x -> x.setReadNegativeStrandFlag(true));

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readStart1, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals(readBases1, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigar1, readInfo.ConsensusRead.getCigarString());

    }

    @Test
    public void testSbxIndelConsensusWithSoftClipsJitter()
    {
        // requires ability to have different unclipped and/or aligned start locations
        List<SAMRecord> reads = Lists.newArrayList();

        // 1:   ACGT AAC CGG TT AC
        //      SSSS MMM MMMDMM MS

        // 2:     GT AACCCGG TT ACGT
        //        SS MMMIMMMDMM MSSS

        // 3:    CGT AACCCGGGTT A
        //       SSM MMMIMMMMMM S

        // net: ACGT AACCCGG TT ACGT
        //      SSSS MMMIMMMDMM MSSS
        // cigar: 3S 3M 1I 3M 1D 3M 3S

        String readBases1 = "ACGTAACCGGTTAC";
        String readCigar1 = "4S6M1D3M1S";
        int readStart1 = 5; // unclipped start is 1

        String readBases2 = "GTAACCCGGTTACGT";
        String readCigar2 = "2S3M1I3M1D3M3S";
        int readStart2 = 5; // unclipped start is 3

        String readBases3 = "CGTAACCCGGGTTA";
        String readCigar3 = "2S4M1I6M1S";
        int readStart3 = 4; // unclipped start is 2

        String readBasesNet = "ACGTAACCCGGTTACGT";
        String readCigarNet = "4S3M1I3M1D3M3S";

        byte[] readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read1 = createSamRecord(readBases1, readStart1, readBaseQuals1, readCigar1);
        reads.add(read1);

        byte[] readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read2 = createSamRecord(readBases2, readStart2, readBaseQuals2, readCigar2);
        reads.add(read2);

        byte[] readBaseQuals3 = buildBaseQuals(readBases3.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read3 = createSamRecord(readBases3, readStart3, readBaseQuals3, readCigar3);
        reads.add(read3);

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readStart1, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals(readBasesNet, readInfo.ConsensusRead.getReadString());
        assertEquals(readCigarNet, readInfo.ConsensusRead.getCigarString());

        reads.clear();

        // on the reverse strand
        // 1:     C AACC C GGTT AAA
        //        S MMMM I MMMM SSS

        // 2:    AC AACC C GGTT TA
        //       SS MMMM I MMMM MS

        // 3:   AAA AACC C GGTT T
        //      SSM MMMM I MMMM S

        // net: AAC AACC C GGTT TAA
        //      SSS MMMM I MMMM MSS

        readBases1 = "CAACCCGGTTAAA";
        readCigar1 = "1S4M1I4M3S";
        readStart1 = 5; // unclipped start is 4

        readBases2 = "ACAACCCGGTTTA";
        readCigar2 = "2S4M1I5M1S";
        readStart2 = 5; // unclipped start is 3

        readBases3 = "AAAAACCCGGTTT";
        readCigar3 = "2S5M1I4M1S";
        readStart3 = 4; // unclipped start is 2

        readBasesNet = "AACAACCCGGTTTAA";
        readCigarNet = "3S4M1I4M3S";

        readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_SIMPLEX_QUAL);
        read1 = createSamRecord(readBases1, readStart1, readBaseQuals1, readCigar1);
        reads.add(read1);

        readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_SIMPLEX_QUAL);
        read2 = createSamRecord(readBases2, readStart2, readBaseQuals2, readCigar2);
        reads.add(read2);

        readBaseQuals3 = buildBaseQuals(readBases3.length(), RAW_SIMPLEX_QUAL);
        read3 = createSamRecord(readBases3, readStart3, readBaseQuals3, readCigar3);
        reads.add(read3);

        reads.forEach(x -> x.setReadNegativeStrandFlag(true));

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readStart2, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals(readCigarNet, readInfo.ConsensusRead.getCigarString());
        assertEquals(readBasesNet, readInfo.ConsensusRead.getReadString());
    }

    @Test
    public void testSbxIndelConsensusWithSoftClipsJitter2()
    {
        List<SAMRecord> reads = Lists.newArrayList();

        // 1:  AACC   GGTT
        //     SSSS   MMMM

        // 2:    CC T GGTT
        //       MM I MMMM

        // 3:     C T GGTT
        //        M I MMMM

        // net: AACC T GGTT
        //      SSMM I MMMM
        // cigar: 2S 2M 1I 4M

        String readBases1 = "AACCGGTT";
        String readCigar1 = "4S4M";
        int readStart1 = 5; // unclipped start is 1

        // pos:              34 5678
        // index:            0123456
        String readBases2 = "CCTGGTT";
        String readCigar2 = "2M1I4M";
        int readStart2 = 3; // unclipped start is 3

        String readBases3 = "CTGGTT";
        String readCigar3 = "1M1I4M";
        int readStart3 = 4; // unclipped start is 4

        String readBasesNet = "AACCTGGTT";
        String readCigarNet = "2S2M1I4M";

        byte[] readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read1 = createSamRecord(readBases1, readStart1, readBaseQuals1, readCigar1);
        read1.setReadNegativeStrandFlag(true);
        reads.add(read1);

        byte[] readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read2 = createSamRecord(readBases2, readStart2, readBaseQuals2, readCigar2);
        read2.setReadNegativeStrandFlag(true);
        reads.add(read2);

        byte[] readBaseQuals3 = buildBaseQuals(readBases3.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read3 = createSamRecord(readBases3, readStart3, readBaseQuals3, readCigar3);
        reads.add(read3);
        reads.forEach(x -> x.setReadNegativeStrandFlag(true));

        ConsensusReadInfo readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readStart2, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals(readCigarNet, readInfo.ConsensusRead.getCigarString());
        assertEquals(readBasesNet, readInfo.ConsensusRead.getReadString());

        // test a delete overlapping soft-clips and prior to aligned read starts
        reads.clear();

        // 1:      T AACC
        //         S MMMM

        // 2:     TT AACC
        //        SS MMMM

        // 3:        AACC
        //           MMMM

        // 4:    GT  AACC
        //       MMD MMMM

        // new:  GT  AACC
        //       MMD MMMM
        // cigar: 2M 1D 4M

        readBases1 = "TAACC";
        readCigar1 = "1S4M";
        readStart1 = 10; // unclipped start is 9

        readBases2 = "TTAACC";
        readCigar2 = "2S4M";
        readStart2 = 10; // unclipped start is 8

        readBases3 = "AACC";
        readCigar3 = "4M";
        readStart3 = 10; // unclipped start is 10

        String readBases4 = "GTAACC";
        String readCigar4 = "2M1D4M";
        int readStart4 = 7; // unclipped start is 7

        readBasesNet = readBases4;
        readCigarNet = readCigar4;

        readBaseQuals1 = buildBaseQuals(readBases1.length(), RAW_SIMPLEX_QUAL);
        read1 = createSamRecord(readBases1, readStart1, readBaseQuals1, readCigar1);
        reads.add(read1);

        readBaseQuals2 = buildBaseQuals(readBases2.length(), RAW_SIMPLEX_QUAL);
        read2 = createSamRecord(readBases2, readStart2, readBaseQuals2, readCigar2);
        reads.add(read2);

        readBaseQuals3 = buildBaseQuals(readBases3.length(), RAW_SIMPLEX_QUAL);
        read3 = createSamRecord(readBases3, readStart3, readBaseQuals3, readCigar3);
        reads.add(read3);

        byte[] readBaseQuals4 = buildBaseQuals(readBases4.length(), RAW_SIMPLEX_QUAL);
        SAMRecord read4 = createSamRecord(readBases4, readStart4, readBaseQuals4, readCigar4);
        reads.add(read4);

        readInfo = createConsensusRead(mConsensusReads, reads, "");
        assertEquals(INDEL_MISMATCH, readInfo.Outcome);
        assertEquals(readStart4, readInfo.ConsensusRead.getAlignmentStart());
        assertEquals(readCigarNet, readInfo.ConsensusRead.getCigarString());
        assertEquals(readBasesNet, readInfo.ConsensusRead.getReadString());
    }

    @Test
    public void testSbxFinaliseReads()
    {
        int position = 1;
        String refBases = REF_BASES.substring(position, 11);

        String readBases = "T" + refBases.substring(1);

        byte[] baseQuals = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);
        SAMRecord read = createSamRecord(readBases, position, baseQuals);

        SbxRoutines.finaliseRead(mRefGenome, read);

        for(int i = 0; i < read.getBaseQualities().length; ++i)
        {
            assertEquals(SBX_DUPLEX_QUAL, read.getBaseQualities()[i]);
        }

        ConsensusType consensusType = extractConsensusType(read);
        int duplexBaseIndex = SbxRoutines.findMaxDuplexBaseIndex(List.of(read));
        assertEquals(ConsensusType.NONE, consensusType);
        assertEquals(0, duplexBaseIndex);

        baseQuals = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);
        read = createSamRecord(readBases, position, baseQuals);
        read.setReadNegativeStrandFlag(true);
        read.setAttribute(CONSENSUS_TYPE_ATTRIBUTE, ConsensusType.SINGLE.toString());
        setDupluxBaseIndex(read);

        SbxRoutines.finaliseRead(mRefGenome, read);

        consensusType = extractConsensusType(read);
        duplexBaseIndex = extractDuplexBaseIndex(read);
        assertEquals(ConsensusType.DUAL, consensusType);
        assertEquals(9, duplexBaseIndex);

        // mark the transition point mid-way through the read
        baseQuals = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);

        int lastIndex = readBases.length() - 1;
        baseQuals[lastIndex] = RAW_SIMPLEX_QUAL;
        baseQuals[lastIndex - 1] = RAW_SIMPLEX_QUAL;
        baseQuals[lastIndex - 2] = RAW_SIMPLEX_QUAL;

        read = createSamRecord(readBases, position, baseQuals);
        read.setAttribute(CONSENSUS_TYPE_ATTRIBUTE, ConsensusType.SINGLE.toString());
        read.setReadNegativeStrandFlag(true);
        setDupluxBaseIndex(read);

        SbxRoutines.finaliseRead(mRefGenome, read);

        consensusType = extractConsensusType(read);
        duplexBaseIndex = extractDuplexBaseIndex(read);
        assertEquals(ConsensusType.DUAL, consensusType);
        assertEquals(6, duplexBaseIndex);

        // test duplex mismatch bases
        baseQuals = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);

        baseQuals[3] = SBX_DUPLEX_MISMATCH_QUAL;
        baseQuals[7] = SBX_DUPLEX_MISMATCH_QUAL;
        baseQuals[8] = SBX_DUPLEX_MISMATCH_QUAL;

        //                x   xx
        //             0123456789
        // ref bases:  AAAAACCCCC
        // read bases: AATAACCCCG
        readBases = "AATAACCCCG";
        read = createSamRecord(readBases, position, baseQuals);
        SbxRoutines.finaliseRead(mRefGenome, read);

        assertEquals(SBX_DUPLEX_ADJACENT_2_3_QUAL, read.getBaseQualities()[0]);
        assertEquals(SBX_DUPLEX_ADJACENT_2_3_QUAL, read.getBaseQualities()[1]);
        assertEquals(SBX_DUPLEX_ADJACENT_1_QUAL, read.getBaseQualities()[2]);
        assertEquals(SBX_DUPLEX_MISMATCH_QUAL, read.getBaseQualities()[3]);
        assertEquals(SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH, read.getBaseQualities()[4]);
        assertEquals(SBX_DUPLEX_ADJACENT_2_3_QUAL, read.getBaseQualities()[5]);
        assertEquals(SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH, read.getBaseQualities()[6]);
        assertEquals(SBX_DUPLEX_MISMATCH_QUAL, read.getBaseQualities()[7]);
        assertEquals(SBX_DUPLEX_MISMATCH_QUAL, read.getBaseQualities()[8]);
        assertEquals(SBX_DUPLEX_ADJACENT_1_QUAL, read.getBaseQualities()[9]);

        // test a conversion of an X to M - not sure if this will occur in actual BAMs
        read = createSamRecord(readBases, position, baseQuals);
        read.setCigarString("2X3I2X1D3X");

        SbxRoutines.finaliseRead(mRefGenome, read);
        assertEquals("2M3I2M1D3M", read.getCigarString());
    }

    @Test
    public void testReplaceXWithMCigar()
    {
        // test a conversion of an X to M - not sure if this will occur in actual BAMs
        int position = 1;
        String readBases = REF_BASES.substring(position, 11);
        byte[] baseQuals = buildBaseQuals(readBases.length(), RAW_DUPLEX_QUAL);

        SAMRecord read = createSamRecord(readBases, position, baseQuals);
        read.setCigarString("10X");

        SbxRoutines.finaliseRead(mRefGenome, read);
        assertEquals("10M", read.getCigarString());

        // multiple
        read = createSamRecord(readBases, position, baseQuals);
        read.setCigarString("2X3I2X1D3X");

        SbxRoutines.finaliseRead(mRefGenome, read);
        assertEquals("2M3I2M1D3M", read.getCigarString());

        // requires collapsing
        read = createSamRecord(readBases, position, baseQuals);
        read.setCigarString("2X2M2X2M2X");

        SbxRoutines.finaliseRead(mRefGenome, read);
        assertEquals("10M", read.getCigarString());
    }

    private static void setDupluxBaseIndex(final SAMRecord record)
    {
        int firstDuplexBaseIndex = SbxRoutines.findMaxDuplexBaseIndex(List.of(record));

        if(firstDuplexBaseIndex >= 0)
            record.setAttribute(SBX_DUPLEX_READ_INDEX_TAG, firstDuplexBaseIndex);
    }

    private static SAMRecord createSamRecord(final String readBases, final int position, final byte[] baseQualities)
    {
        String cigar = format("%dM", readBases.length());
        return createSamRecord(readBases, position, baseQualities, cigar);
    }

    private static SAMRecord createSamRecord(final String readBases, final int position, final byte[] baseQualities, final String cigar)
    {
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, position, readBases, cigar, CHR_1, 300,
                false, false, null, true, TEST_READ_CIGAR);
        record.setBaseQualities(baseQualities);
        return record;
    }
}
