package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsToStr;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractConsensusType;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.RAW_SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_ADJACENT_1_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_ADJACENT_2_3_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_READ_INDEX_TAG;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.extractDuplexBaseIndex;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildBaseQuals;
import static com.hartwig.hmftools.redux.consensus.SbxConsensusTest.createSamRecord;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.correctInvalidCigar;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.ensureMirroredMismatchRepeatQuals;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.hasInvalidCigar;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.redux.ReduxConfig;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class SbxFinaliseReadTest
{
    private final MockRefGenome mRefGenome;
    private final ConsensusReads mConsensusReads;

    //                                                 20        30        40        50
    //                                        12345678901234567890123456789012345678901234567890
    private static final String REF_BASES = "XAAGAGCTCCAGGGGGTTTTTACGTACGTACGTAACCGGTTAACCGGTT";

    public SbxFinaliseReadTest()
    {
        mRefGenome = new MockRefGenome(true);
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
        // ref bases:  AAGAGCTCCA
        // read bases: AATAGCTCCG
        readBases =   "AATAGCTCCG";
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

    @Test
    public void testMismatchRepeatMirroring()
    {
        //                           10        20        30
        //                 0123456789012345678901234567890123456789
        String basesStr = "TTTTAAAATTTACACACTTTAAGAAGTTT";
        byte[] baseQuals = buildBaseQuals(basesStr.length(), SBX_DUPLEX_QUAL);

        // test left to right repeats
        markLowQualBases(baseQuals, List.of(4, 11, 12, 20, 21, 22));

        ensureMirroredMismatchRepeatQuals(0, basesStr.length(), basesStr.getBytes(), baseQuals);

        checkBases(baseQuals, List.of(4, 7, 11, 12, 15, 16, 20, 21, 22, 23, 24, 25));

        baseQuals = buildBaseQuals(basesStr.length(), SBX_DUPLEX_QUAL);
        markLowQualBases(baseQuals, List.of(7, 15, 16, 23, 24, 25));

        ensureMirroredMismatchRepeatQuals(0, basesStr.length(), basesStr.getBytes(), baseQuals);

        checkBases(baseQuals, List.of(4, 7, 11, 12, 15, 16, 20, 21, 22, 23, 24, 25));
    }

    @Test
    public void testFixCigars()
    {
        List<CigarElement> cigarElements = CigarUtils.cigarElementsFromStr("10I20M");

        assertTrue(hasInvalidCigar(cigarElements));

        correctInvalidCigar(cigarElements);

        assertEquals("10S20M", cigarElementsToStr(cigarElements));

        cigarElements = CigarUtils.cigarElementsFromStr("20M10I");

        assertTrue(hasInvalidCigar(cigarElements));

        correctInvalidCigar(cigarElements);

        assertEquals("20M10S", cigarElementsToStr(cigarElements));

        cigarElements = CigarUtils.cigarElementsFromStr("10I20M10I");

        assertTrue(hasInvalidCigar(cigarElements));

        correctInvalidCigar(cigarElements);

        assertEquals("10S20M10S", cigarElementsToStr(cigarElements));

        cigarElements = CigarUtils.cigarElementsFromStr("10S10I20M10I10S");

        assertTrue(hasInvalidCigar(cigarElements));

        correctInvalidCigar(cigarElements);

        assertEquals("20S20M20S", cigarElementsToStr(cigarElements));
    }

    private static void markLowQualBases(final byte[] baseQuals, final List<Integer> lowQualBases)
    {
        for(int i = 0; i < baseQuals.length; ++i)
        {
            baseQuals[i] = SBX_DUPLEX_MISMATCH_QUAL;
        }
    }

    private static boolean checkBases(final byte[] baseQuals, final List<Integer> lowQualBases)
    {
        for(int i = 0; i < baseQuals.length; ++i)
        {
            if(lowQualBases.contains(i))
            {
                if(baseQuals[i] != SBX_DUPLEX_MISMATCH_QUAL)
                    return false;
            }
            else
            {
                if(baseQuals[i] == SBX_DUPLEX_MISMATCH_QUAL)
                    return false;
            }
        }

        return true;
    }
}
