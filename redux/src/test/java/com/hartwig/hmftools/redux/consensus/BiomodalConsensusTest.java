package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.BASE_MODIFICATIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.BASE_QUAL_MINIMUM;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.MM_PREFIX;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.MM_SUFFIX;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.getMMValueFromModCReadIndices;
import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.redux.consensus.ConsensusOutcome.INDEL_MISMATCH;

import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.consensus.ConsensusOutcome;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;
import com.hartwig.hmftools.redux.consensus.RefGenome;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class BiomodalConsensusTest
{
    private static final char MODC_BASE = 'c';
    private static final String QUAL_25 = String.valueOf(phredToFastq(25));

    private final RefGenome mRefGenome;
    private final ConsensusReads mConsensusReads;
    private final Map<Byte, Byte> mNextBaseMap;

    public BiomodalConsensusTest()
    {
        MockRefGenome mockRefGenome = new MockRefGenome(true);
        mockRefGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        mockRefGenome.ChromosomeLengths.put(CHR_1, REF_BASES.length());

        mConsensusReads = new ConsensusReads(mockRefGenome, BIOMODAL);
        mRefGenome = new RefGenome(mockRefGenome);

        mNextBaseMap = Maps.newHashMap();
        mNextBaseMap.put((byte) 'G', (byte) 'C');
        mNextBaseMap.put((byte) 'C', (byte) 'A');
        mNextBaseMap.put((byte) 'A', (byte) 'T');
        mNextBaseMap.put((byte) 'T', (byte) 'G');
    }

    /*
    @Test
    public void testBiomodalConsensusReadSimple()
    {
        int alignmentStart = 20;
        int readLength = 50;
        String readStr = REF_BASES.substring(alignmentStart - 1, alignmentStart - 1 + readLength);
        String qualStr = QUAL_25.repeat(readLength);
        String cigar = readLength + "M";

        SAMRecord read1 = createBiomodalSamRecord("READ_001", CHR_1, alignmentStart, readStr, qualStr, cigar, true);
        SAMRecord read2 = createBiomodalSamRecord("READ_002", CHR_1, alignmentStart, readStr, qualStr, cigar, true);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);

        ConsensusReadInfo consensusOutput = mConsensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        assertEquals(ALIGNMENT_ONLY, consensusOutcome);
        assertEquals(alignmentStart, consensusRead.getAlignmentStart());
        assertEquals(cigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertEquals(MM_PREFIX + MM_SUFFIX, consensusRead.getStringAttribute(BASE_MODIFICATIONS_ATTRIBUTE));
    }

    @Test
    public void testBiomodalConsensusReadWithClippingAndIndels()
    {
        // construct read1 and read2
        String qualStr = QUAL_25.repeat(50);
        StringBuilder readStr = new StringBuilder();
        String cigar = "1S10M1D10M1I28M";

        readStr.append((char) ((byte) mNextBaseMap.get(mRefGenome.getRefBase(CHR_1, 20))));
        for(int i = 21; i < 31; i++)
            readStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        for(int i = 32; i < 42; i++)
            readStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        readStr.append((char) mRefGenome.getRefBase(CHR_1, 42));
        for(int i = 42; i < 70; i++)
            readStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        assertEquals(qualStr.length(), readStr.length());
        SAMRecord read1 = createBiomodalSamRecord("READ_001", CHR_1, 21, readStr.toString(), qualStr, cigar, true);
        SAMRecord read2 = createBiomodalSamRecord("READ_002", CHR_1, 21, readStr.toString(), qualStr, cigar, true);

        // construct read3
        readStr = new StringBuilder();
        cigar = "10M1D40M";

        for(int i = 21; i < 31; i++)
            readStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        for(int i = 32; i < 72; i++)
            readStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        assertEquals(qualStr.length(), readStr.length());
        SAMRecord read3 = createBiomodalSamRecord("READ_003", CHR_1, 21, readStr.toString(), qualStr, cigar, true);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = mConsensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        int expectedAlignmentStart = 21;
        String expectedCigar = "1S10M1D10M1I30M";
        StringBuilder expectedReadStr = new StringBuilder();
        expectedReadStr.append((char) ((byte) mNextBaseMap.get(mRefGenome.getRefBase(CHR_1, 20))));
        for(int i = 21; i < 31; i++)
            expectedReadStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        for(int i = 32; i < 42; i++)
            expectedReadStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        expectedReadStr.append((char) mRefGenome.getRefBase(CHR_1, 42));
        for(int i = 42; i < 72; i++)
            expectedReadStr.append((char) mRefGenome.getRefBase(CHR_1, i));

        String expectedQualStr = QUAL_25.repeat(expectedReadStr.length());

        assertEquals(INDEL_MISMATCH, consensusOutcome);
        assertEquals(expectedAlignmentStart, consensusRead.getAlignmentStart());
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr.toString(), consensusRead.getReadString());
        assertEquals(expectedQualStr, consensusRead.getBaseQualityString());
        assertEquals(MM_PREFIX + MM_SUFFIX, consensusRead.getStringAttribute(BASE_MODIFICATIONS_ATTRIBUTE));
    }

    @Test
    public void testBiomodalConsensusModCForward()
    {
        int alignmentStart = 20;
        int readLength = 50;
        String readStr = REF_BASES.substring(alignmentStart - 1, alignmentStart - 1 + readLength);
        String qualStr = QUAL_25.repeat(readLength);
        String cigar = readLength + "M";

        StringBuilder[] modCReadStrs = new StringBuilder[] { new StringBuilder(readStr), new StringBuilder(readStr), new StringBuilder(readStr) };
        List<Integer> cReadIndices = Lists.newArrayList();
        for(int i = 0; i < readStr.length(); i++)
        {
            if(readStr.charAt(i) == 'C')
                cReadIndices.add(i);
        }

        assertTrue(cReadIndices.size() >= 4);

        SortedSet<Integer> consensusModCReadIndices = Sets.newTreeSet();
        for(int i = 0; i < cReadIndices.size(); i++)
        {
            int readIndex = cReadIndices.get(i);
            int modCCount = i % 4;
            for(int j = 0; j < modCCount; j++)
                modCReadStrs[j].setCharAt(readIndex, MODC_BASE);

            if(modCCount >= 2)
                consensusModCReadIndices.add(readIndex);
        }

        SAMRecord read1 = createBiomodalSamRecord("READ_001", CHR_1, alignmentStart, modCReadStrs[0].toString(), qualStr, cigar, true);
        SAMRecord read2 = createBiomodalSamRecord("READ_002", CHR_1, alignmentStart, modCReadStrs[1].toString(), qualStr, cigar, true);
        SAMRecord read3 = createBiomodalSamRecord("READ_003", CHR_1, alignmentStart, modCReadStrs[2].toString(), qualStr, cigar, true);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);

        ConsensusReadInfo consensusOutput = mConsensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;
        String expectedMMValue = getMMValueFromModCReadIndices(consensusRead.getReadBases(), consensusModCReadIndices, true);

        assertEquals(ALIGNMENT_ONLY, consensusOutcome);
        assertEquals(alignmentStart, consensusRead.getAlignmentStart());
        assertEquals(cigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertFalse(consensusRead.getReadNegativeStrandFlag());
        assertEquals(expectedMMValue, consensusRead.getStringAttribute(BASE_MODIFICATIONS_ATTRIBUTE));
    }

    @Test
    public void testBiomodalConsensusModCReverse()
    {
        int alignmentStart = 20;
        int readLength = 50;
        String readStr = REF_BASES.substring(alignmentStart - 1, alignmentStart - 1 + readLength);
        String qualStr = QUAL_25.repeat(readLength);
        String cigar = readLength + "M";

        StringBuilder[] modCReadStrs = new StringBuilder[] { new StringBuilder(readStr), new StringBuilder(readStr), new StringBuilder(readStr) };
        List<Integer> gReadIndices = Lists.newArrayList();
        for(int i = 0; i < readStr.length(); i++)
        {
            if(readStr.charAt(i) == 'G')
                gReadIndices.add(i);
        }

        assertTrue(gReadIndices.size() >= 4);

        SortedSet<Integer> consensusModCReadIndices = Sets.newTreeSet();
        for(int i = 0; i < gReadIndices.size(); i++)
        {
            int readIndex = gReadIndices.get(i);
            int modCCount = i % 4;
            for(int j = 0; j < modCCount; j++)
                modCReadStrs[j].setCharAt(readIndex, MODC_BASE);

            if(modCCount >= 2)
                consensusModCReadIndices.add(readIndex);
        }

        SAMRecord read1 = createBiomodalSamRecord("READ_001", CHR_1, alignmentStart, modCReadStrs[0].toString(), qualStr, cigar, false);
        SAMRecord read2 = createBiomodalSamRecord("READ_002", CHR_1, alignmentStart, modCReadStrs[1].toString(), qualStr, cigar, false);
        SAMRecord read3 = createBiomodalSamRecord("READ_003", CHR_1, alignmentStart, modCReadStrs[2].toString(), qualStr, cigar, false);
        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);

        ConsensusReadInfo consensusOutput = mConsensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;
        String expectedMMValue = getMMValueFromModCReadIndices(consensusRead.getReadBases(), consensusModCReadIndices, false);

        assertEquals(ALIGNMENT_ONLY, consensusOutcome);
        assertEquals(alignmentStart, consensusRead.getAlignmentStart());
        assertEquals(cigar, consensusRead.getCigarString());
        assertEquals(readStr, consensusRead.getReadString());
        assertEquals(qualStr, consensusRead.getBaseQualityString());
        assertTrue(consensusRead.getReadNegativeStrandFlag());
        assertEquals(expectedMMValue, consensusRead.getStringAttribute(BASE_MODIFICATIONS_ATTRIBUTE));
    }

    @Test
    public void testBiomodalConsensusDropNonStrictMajorityInserts()
    {
        String refBases = "A".repeat(20) + "AAAAAAAA" + "A".repeat(20);
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr1 = "A".repeat(20) + "AATAATAATAA" + "A".repeat(20);
        String qualStr1 = QUAL_25.repeat(readStr1.length());
        String cigar1 = "22M1I2M1I2M1I22M";

        String readStr2 = "A".repeat(20) + "AAAATAATAA" + "A".repeat(20);
        String qualStr2 = QUAL_25.repeat(readStr2.length());
        String cigar2 = "24M1I2M1I22M";

        String readStr3 = "A".repeat(20) + "AAAAAATAA" + "A".repeat(20);
        String qualStr3 = QUAL_25.repeat(readStr3.length());
        String cigar3 = "26M1I22M";

        SAMRecord read1 = createBiomodalSamRecord("READ_001", CHR_1, 1, readStr1, qualStr1, cigar1, true);
        SAMRecord read2 = createBiomodalSamRecord("READ_002", CHR_1, 1, readStr2, qualStr2, cigar2, true);
        SAMRecord read3 = createBiomodalSamRecord("READ_003", CHR_1, 1, readStr3, qualStr3, cigar3, true);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        int expectedAlignmentStart = 1;
        String expectedCigar = cigar2;
        String expectedReadStr = readStr2;
        String expectedQualStr = qualStr2;

        assertEquals(INDEL_MISMATCH, consensusOutcome);
        assertEquals(expectedAlignmentStart, consensusRead.getAlignmentStart());
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr, consensusRead.getReadString());
        assertEquals(expectedQualStr, consensusRead.getBaseQualityString());
    }

    @Test
    public void testBiomodalConsensusKeepStrictMajorityDels()
    {
        String refBases = "A".repeat(20) + "ATGCATGC" + "A".repeat(20);
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr1 = "A".repeat(20) + "AGCTC" + "A".repeat(20);
        String qualStr1 = QUAL_25.repeat(readStr1.length());
        String cigar1 = "21M1D2M1D1M1D21M";

        String readStr2 = "A".repeat(20) + "ATGCTC" + "A".repeat(20);
        String qualStr2 = QUAL_25.repeat(readStr2.length());
        String cigar2 = "24M1D1M1D21M";

        String readStr3 = "A".repeat(20) + "ATGCATC" + "A".repeat(20);
        String qualStr3 = QUAL_25.repeat(readStr3.length());
        String cigar3 = "26M1D21M";

        SAMRecord read1 = createBiomodalSamRecord("READ_001", CHR_1, 1, readStr1, qualStr1, cigar1, true);
        SAMRecord read2 = createBiomodalSamRecord("READ_002", CHR_1, 1, readStr2, qualStr2, cigar2, true);
        SAMRecord read3 = createBiomodalSamRecord("READ_003", CHR_1, 1, readStr3, qualStr3, cigar3, true);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2, read3);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        int expectedAlignmentStart = 1;
        String expectedCigar = cigar2;
        String expectedReadStr = readStr2;
        String expectedQualStr = qualStr2;

        assertEquals(INDEL_MISMATCH, consensusOutcome);
        assertEquals(expectedAlignmentStart, consensusRead.getAlignmentStart());
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr, consensusRead.getReadString());
        assertEquals(expectedQualStr, consensusRead.getBaseQualityString());
    }

    @Test
    public void testBiomodalConsensusNoReplacementOfSoftClipWithRef()
    {
        String refBases = "AA";
        ConsensusReads consensusReads = getConsensusReads(refBases);

        String readStr1 = "TA";
        String readStr2 = "CA";

        String qualStr = QUAL_25.repeat(readStr1.length());
        String cigar = "1S1M";

        SAMRecord read1 = createBiomodalSamRecord("READ_001", CHR_1, 2, readStr1, qualStr, cigar, true);
        SAMRecord read2 = createBiomodalSamRecord("READ_002", CHR_1, 2, readStr2, qualStr, cigar, true);

        List<SAMRecord> reads = Lists.newArrayList(read1, read2);
        FragmentCoords coords = FragmentCoords.fromRead(read1, false);
        ConsensusReadInfo consensusOutput = consensusReads.createConsensusRead(reads, coords, null);
        ConsensusOutcome consensusOutcome = consensusOutput.Outcome;
        SAMRecord consensusRead = consensusOutput.ConsensusRead;

        int expectedAlignmentStart = 2;
        String expectedCigar = "1S1M";
        String expectedReadStr = "TA";
        String expectedQualStr = phredToFastq(BASE_QUAL_MINIMUM) + QUAL_25;

        assertEquals(ALIGNMENT_ONLY, consensusOutcome);
        assertEquals(expectedAlignmentStart, consensusRead.getAlignmentStart());
        assertEquals(expectedCigar, consensusRead.getCigarString());
        assertEquals(expectedReadStr, consensusRead.getReadString());
        assertEquals(expectedQualStr, consensusRead.getBaseQualityString());
    }
    */

    private static SAMRecord createBiomodalSamRecord(
            final String readName, final String chromosome, int alignmentStart, final String modCReadStr, final String qualStr,
            final String cigar, boolean isForward)
    {
        char modCBase = isForward ? 'C' : swapDnaBase('C');
        StringBuilder readStr = new StringBuilder(modCReadStr);
        SortedSet<Integer> modCReadIndices = Sets.newTreeSet();
        for(int i = 0; i < readStr.length(); i++)
        {
            char base = readStr.charAt(i);
            if(base != MODC_BASE)
                continue;

            readStr.setCharAt(i, modCBase);
            modCReadIndices.add(i);
        }

        String mmValue = getMMValueFromModCReadIndices(readStr.toString().getBytes(), modCReadIndices, isForward);

        SAMRecord read = createSamRecordUnpaired(readName, chromosome, alignmentStart, readStr.toString(), cigar, !isForward, false, null);
        read.setBaseQualityString(qualStr);
        read.setMappingQuality(60);
        read.setAttribute(BASE_MODIFICATIONS_ATTRIBUTE, mmValue);

        return read;
    }

    private static ConsensusReads getConsensusReads(final String refBases)
    {
        MockRefGenome mockRefGenome = new MockRefGenome(true);
        mockRefGenome.RefGenomeMap.put(CHR_1, refBases);
        mockRefGenome.ChromosomeLengths.put(CHR_1, refBases.length());
        return new ConsensusReads(mockRefGenome, BIOMODAL);
    }
}
