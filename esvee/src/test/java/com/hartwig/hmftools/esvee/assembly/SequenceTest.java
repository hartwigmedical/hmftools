package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.buildFlags;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.setSbxSequencing;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.buildTrimmedRefBaseSequence;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findDualBaseRepeat;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findDualDualRepeat;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findMultiBaseRepeat;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findRepeats;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findSingleBaseRepeat;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findSingleOrDualRepeat;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findTripleBaseRepeat;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.assembly.read.Read;

import org.junit.Test;

public class SequenceTest
{
    @Test
    public void testSeqBuildSnvMismatches()
    {
        String refBuffer = "CCCCC";
        //                      012345678901
        String readExtBases1 = "ACGTAAGTACGT";
        String readBases1 = refBuffer + readExtBases1;
        String readBases2 = readBases1;
        String readExtBases3 = "ACGGAAGGACGG";
        String readBases3 = refBuffer + readExtBases3;

        String cigar = makeCigarString(readBases1, 0, 0);

        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases1, cigar);
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases2, cigar);
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases3, cigar);

        int consensusBaseLength = readBases1.length() - refBuffer.length();
        boolean buildForwards = true;

        List<ReadParseState> readParseStates = Lists.newArrayList(
                new ReadParseState(buildForwards, read1, refBuffer.length()),
                new ReadParseState(buildForwards, read2, refBuffer.length()),
                new ReadParseState(buildForwards, read3, refBuffer.length()));

        SequenceBuilder seqBuilder = new SequenceBuilder(readParseStates, buildForwards, consensusBaseLength);

        String consensusStr = seqBuilder.baseString();
        assertEquals(readExtBases1, consensusStr);
        assertTrue(seqBuilder.repeats().isEmpty());

        ReadParseState readState3 = readParseStates.get(2);
        assertEquals(3, readState3.mismatchInfos().size());

        // test favouring high over medium quals
        setSbxSequencing();

        read1.getBaseQuality()[8] = SbxBamUtils.SBX_SIMPLEX_QUAL;
        read2.getBaseQuality()[8] = SbxBamUtils.SBX_DUPLEX_QUAL;
        read3.getBaseQuality()[8] = SbxBamUtils.SBX_DUPLEX_QUAL;

        read1.getBaseQuality()[12] = SbxBamUtils.SBX_SIMPLEX_QUAL;
        read2.getBaseQuality()[12] = SbxBamUtils.SBX_SIMPLEX_QUAL;
        read3.getBaseQuality()[12] = SbxBamUtils.SBX_DUPLEX_QUAL;

        readParseStates.forEach(x -> x.resetAll());
        seqBuilder = new SequenceBuilder(readParseStates, buildForwards, consensusBaseLength);

        consensusStr = seqBuilder.baseString();
        assertEquals("ACGTAAGGACGT", consensusStr);

        setIlluminaSequencing();

        // repeat test building the other direction
        readBases1 = readExtBases1 + refBuffer;
        readBases2 = readBases1;
        readBases3 = readExtBases3 + refBuffer;

        read1 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases1, cigar);
        read2 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases2, cigar);
        read3 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases3, cigar);

        buildForwards = false;

        readParseStates = Lists.newArrayList(
                new ReadParseState(buildForwards, read1, readExtBases1.length() - 1),
                new ReadParseState(buildForwards, read2, readExtBases1.length() - 1),
                new ReadParseState(buildForwards, read3, readExtBases1.length() - 1));

        seqBuilder = new SequenceBuilder(readParseStates, buildForwards, consensusBaseLength);

        consensusStr = seqBuilder.baseString();
        assertEquals(readExtBases1, consensusStr);
        assertTrue(seqBuilder.repeats().isEmpty());
    }

    @Test
    public void testSeqBuildIndelMismatches()
    {
        String refBuffer = "CCCCC";
        //                      012345678901
        String readExtBases1 = "ACGTAAGTACGTGGCC";
        String readBases1 = refBuffer + readExtBases1;
        String readBases2 = readBases1;

        // with a del then an insert
        String readExtBases3 = readExtBases1.substring(0, 4) + readExtBases1.substring(5, 10) + "A" + readExtBases1.substring(10);
        String readBases3 = refBuffer + readExtBases3;

        // with a 2 inserts then a del
        String readExtBases4 = readExtBases1.substring(0, 3) + "A" + readExtBases1.substring(3, 6) + "A"
                + readExtBases1.substring(6, 10) + readExtBases1.substring(11);
        String readBases4 = refBuffer + readExtBases4;

        String cigar = makeCigarString(readBases1, 0, 0);

        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases1, cigar);
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases2, cigar);
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases3, cigar);
        Read read4 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases4, cigar);

        boolean buildForwards = true;

        List<ReadParseState> readParseStates = Lists.newArrayList(
                new ReadParseState(buildForwards, read1, refBuffer.length()),
                new ReadParseState(buildForwards, read2, refBuffer.length()),
                new ReadParseState(buildForwards, read3, refBuffer.length()),
                new ReadParseState(buildForwards, read4, refBuffer.length()));

        int maxReadLength = readParseStates.stream().mapToInt(x -> x.read().basesLength()).max().orElse(0);
        int consensusBaseLength = maxReadLength - refBuffer.length();

        SequenceBuilder seqBuilder = new SequenceBuilder(readParseStates, buildForwards, consensusBaseLength);

        String consensusStr = seqBuilder.baseString();
        assertEquals(readExtBases1, consensusStr);
        assertTrue(seqBuilder.repeats().isEmpty());

        ReadParseState readState3 = readParseStates.get(2);
        assertEquals(2, readState3.mismatchInfos().size());
        assertEquals(SequenceDiffType.DELETE, readState3.mismatchInfos().get(0).Type);
        assertEquals(SequenceDiffType.INSERT, readState3.mismatchInfos().get(1).Type);

        ReadParseState readState4 = readParseStates.get(3);
        assertEquals(3, readState4.mismatchInfos().size());
        assertEquals(SequenceDiffType.INSERT, readState4.mismatchInfos().get(0).Type);
        assertEquals(SequenceDiffType.INSERT, readState4.mismatchInfos().get(1).Type);
        assertEquals(SequenceDiffType.DELETE, readState4.mismatchInfos().get(2).Type);
    }

    @Test
    public void testSeqBuildRepeatMismatches1()
    {
        String refBuffer = "CCCCC";
        //                      012345678901
        String readPreBases = "ACGTAACCT";
        String readPostBases = "TGCAGGGG";

        List<String> repeats = List.of("A", "AT", "AAT", "GGCC", "AATTT");
        List<Integer> readRepeatCounts = List.of(6, 6, 4, 5, 7, 8);

        for(int i = 0; i <= 1; ++i)
        {
            boolean buildForwards = (i == 0);

            for(int test = 0; test < repeats.size(); ++test)
            {
                String repeat = repeats.get(test);

                List<Read> reads = Lists.newArrayList();
                List<ReadParseState> readStates = Lists.newArrayList();

                for(int r = 0; r < readRepeatCounts.size(); ++r)
                {
                    int repeatCount = readRepeatCounts.get(r);

                    String readBases = readPreBases + repeat.repeat(repeatCount) + readPostBases;

                    if(buildForwards)
                        readBases = refBuffer + readBases;
                    else
                        readBases += refBuffer;

                    String cigar = makeCigarString(readBases, 0, 0);
                    Read read = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, cigar);
                    reads.add(read);

                    int readIndexStart = buildForwards ? refBuffer.length() : readBases.length() - refBuffer.length() - 1;
                    readStates.add(new ReadParseState(buildForwards, read, readIndexStart));
                }

                int maxReadLength = readStates.stream().mapToInt(x -> x.read().basesLength()).max().orElse(0);
                int consensusBaseLength = maxReadLength - refBuffer.length();

                SequenceBuilder seqBuilder = new SequenceBuilder(readStates, buildForwards, consensusBaseLength);

                String consensusStr = seqBuilder.baseString();

                String firstReadBases = reads.get(0).getBasesString();
                String consensusBases = buildForwards ?
                        firstReadBases.substring(refBuffer.length()) : firstReadBases.substring(0, firstReadBases.length() - refBuffer.length());

                assertEquals(consensusBases, consensusStr);
                assertEquals(1, seqBuilder.repeats().size());

                /*
                ReadParseState readState3 = readParseStates.get(2);
                assertEquals(2, readState3.mismatchInfos().size());
                assertEquals(SequenceDiffType.DELETE, readState3.mismatchInfos().get(0).Type);
                */

            }
        }
    }

    @Test
    public void testRepeatTypes()
    {
        // single base repeats
        //              0123456789
        String bases = "AAACCTTTTT";

        // first check limits
        RepeatInfo repeatInfo = findSingleBaseRepeat(bases.getBytes(), 5);
        assertNotNull(repeatInfo);
        assertEquals("T", repeatInfo.Bases);
        assertEquals(5, repeatInfo.Count);

        repeatInfo = findSingleBaseRepeat(bases.getBytes(), 6);
        assertNotNull(repeatInfo);
        assertEquals("T", repeatInfo.Bases);
        assertEquals(4, repeatInfo.Count);

        repeatInfo = findSingleBaseRepeat(bases.getBytes(), 7);
        assertNull(repeatInfo);

        //       01234567890123456789
        bases = "AAACCTTTTGAAAAAATGC";

        repeatInfo = findSingleBaseRepeat(bases.getBytes(), 0);
        assertNull(repeatInfo);


        repeatInfo = findSingleBaseRepeat(bases.getBytes(), 5);
        assertNotNull(repeatInfo);
        assertEquals("T", repeatInfo.Bases);
        assertEquals(4, repeatInfo.Count);
        assertEquals(5, repeatInfo.Index);

        repeatInfo = findSingleBaseRepeat(bases.getBytes(), 10);
        assertNotNull(repeatInfo);
        assertEquals("A", repeatInfo.Bases);
        assertEquals(6, repeatInfo.Count);
        assertEquals(10, repeatInfo.Index);

        //       01234567890123456789
        bases = "ATATGATATATGGAT";
        repeatInfo = findDualBaseRepeat(bases.getBytes(), 0);
        assertNull(repeatInfo);

        repeatInfo = findDualBaseRepeat(bases.getBytes(), 5);
        assertNotNull(repeatInfo);
        assertEquals("AT", repeatInfo.Bases);
        assertEquals(3, repeatInfo.Count);
        assertEquals(5, repeatInfo.Index);

        //       012345678901234567890123456789
        bases = "ATGATGGGCCATGATGATGATG";
        repeatInfo = findTripleBaseRepeat(bases.getBytes(), 0);
        repeatInfo = findMultiBaseRepeat(bases.getBytes(), 0, 3);

        assertNotNull(repeatInfo);
        assertEquals("ATG", repeatInfo.Bases);
        assertEquals(2, repeatInfo.Count);
        assertEquals(0, repeatInfo.Index);

        repeatInfo = findTripleBaseRepeat(bases.getBytes(), 10);
        repeatInfo = findMultiBaseRepeat(bases.getBytes(), 10, 3);
        assertNotNull(repeatInfo);
        assertEquals("ATG", repeatInfo.Bases);
        assertEquals(4, repeatInfo.Count);
        assertEquals(10, repeatInfo.Index);

        //       012345678901234567890123456789
        bases = "AATTAAGGCCTTCCTTCCTTCGTT";
        repeatInfo = findDualDualRepeat(bases.getBytes(), 0);
        assertNull(repeatInfo);

        repeatInfo = findDualDualRepeat(bases.getBytes(), 8);
        assertNotNull(repeatInfo);
        assertEquals("CCTT", repeatInfo.Bases);
        assertEquals(3, repeatInfo.Count);
        assertEquals(8, repeatInfo.Index);

        // five-base repeat
        bases = "ACTGTACTGTACTGTACTGTAATTGGTTGGAAGGT";
        repeatInfo = findMultiBaseRepeat(bases.getBytes(), 0, 5);
        assertNotNull(repeatInfo);
        assertEquals("ACTGT", repeatInfo.Bases);
        assertEquals(4, repeatInfo.Count);

        // test them all together
        //       0123456789012345678901234567890123456789
        bases = "ATTTTTAACTCTCTAAAGTCGTCGTCAACCTTCCTTCCTTCGTT";
        List<RepeatInfo> repeats = findRepeats(bases.getBytes());
        assertNotNull(repeats);
        assertEquals(4, repeats.size());

        assertEquals("T", repeats.get(0).Bases);
        assertEquals(5, repeats.get(0).Count);
        assertEquals(1, repeats.get(0).Index);

        assertEquals("CT", repeats.get(1).Bases);
        assertEquals(3, repeats.get(1).Count);
        assertEquals(8, repeats.get(1).Index);

        assertEquals("GTC", repeats.get(2).Bases);
        assertEquals(3, repeats.get(2).Count);
        assertEquals(17, repeats.get(2).Index);

        assertEquals("CCTT", repeats.get(3).Bases);
        assertEquals(3, repeats.get(3).Count);
        assertEquals(28, repeats.get(3).Index);
    }

    @Test
    public void testSingleRepeatSequence()
    {
        //              012345678901234567890123
        String bases = "ATTTTTACTGTGTGTGTCAAAAAT";

        RepeatInfo repeatInfo = findSingleOrDualRepeat(bases.getBytes(), 1, true);

        assertNotNull(repeatInfo);
        assertEquals(5, repeatInfo.Count);
        assertEquals("T", repeatInfo.Bases);

        repeatInfo = findSingleOrDualRepeat(bases.getBytes(), 5, false);

        assertNotNull(repeatInfo);
        assertEquals(5, repeatInfo.Count);
        assertEquals("T", repeatInfo.Bases);

        // too few
        repeatInfo = findSingleOrDualRepeat(bases.getBytes(), 4, true);
        assertNull(repeatInfo);

        repeatInfo = findSingleOrDualRepeat(bases.getBytes(), 20, false);
        assertNull(repeatInfo);

        // at the ends
        repeatInfo = findSingleOrDualRepeat(bases.getBytes(), 21, true);
        assertNull(repeatInfo);

        repeatInfo = findSingleOrDualRepeat(bases.getBytes(), 2, false);
        assertNull(repeatInfo);

        // dual repeats
        repeatInfo = findSingleOrDualRepeat(bases.getBytes(), 8, true);

        assertNotNull(repeatInfo);
        assertEquals("TG", repeatInfo.Bases);
        assertEquals(4, repeatInfo.Count);

        repeatInfo = findSingleOrDualRepeat(bases.getBytes(), 15, false);

        assertNotNull(repeatInfo);
        assertEquals("TG", repeatInfo.Bases);
        assertEquals(4, repeatInfo.Count);
    }

    private static final int NO_MISMATCH_LIMIT = -1;

    @Test
    public void testSequenceRepeatComparisons()
    {
        // test where the repeat multiple is just sufficient to count for one of the sequences only

        //                              10           20
        //                   01234 567890 1234 5678 9012
        String firstBases = "ACGTA CTCTCT ACGT GAGA ACGT";

        String secondBases = "ACGTA CTCT ACGT GAGAGA ACGT";
        //                    01234 5678 9012 345678 9012
        //                                10          20

        firstBases = firstBases.replaceAll(" ", "");
        secondBases = secondBases.replaceAll(" ", "");

        List<RepeatInfo> firstRepeats = findRepeats(firstBases.getBytes());
        assertEquals(1, firstRepeats.size());

        byte[] firstBaseQuals = buildDefaultBaseQuals(firstBases.length());

        List<RepeatInfo> secondRepeats = findRepeats(secondBases.getBytes());
        assertEquals(1, secondRepeats.size());

        byte[] secondBaseQuals = buildDefaultBaseQuals(secondBases.length());

        int mismatches = SequenceCompare.compareSequences(
                firstBases.getBytes(), firstBaseQuals, 0, firstBaseQuals.length - 1, firstRepeats,
                secondBases.getBytes(), secondBaseQuals, 0, secondBaseQuals.length - 1, secondRepeats,
                NO_MISMATCH_LIMIT);

        assertEquals(2, mismatches);

        mismatches = SequenceCompare.compareSequences(
                firstBases.getBytes(), firstBaseQuals, 0, firstBaseQuals.length - 1, firstRepeats,
                secondBases.getBytes(), secondBaseQuals, 0, secondBaseQuals.length - 1, secondRepeats,
                NO_MISMATCH_LIMIT);

        assertEquals(2, mismatches);
    }

        @Test
    public void testSequenceComparisons()
    {
        //                                   10          20            30
        //                    012345     6789012345 678 90 12345678   9012
        String firstBases =  "ATTTTT     AACTCTCTCT AAA GG CTGACGTA   TTCC";

        String secondBases = "ATTTTTTT   AACTCTCT   AAA    CTGACGTA G TTCC";
        //                    0123456789 01234567   890    12345678 9 0123
        //                               10           20              30

        firstBases = firstBases.replaceAll(" ", "");
        secondBases = secondBases.replaceAll(" ", "");

        List<RepeatInfo> firstRepeats = findRepeats(firstBases.getBytes());
        assertEquals(2, firstRepeats.size());

        byte[] firstBaseQuals = buildDefaultBaseQuals(firstBases.length());

        List<RepeatInfo> secondRepeats = findRepeats(secondBases.getBytes());
        assertEquals(2, secondRepeats.size());

        byte[] secondBaseQuals = buildDefaultBaseQuals(secondBases.length());

        // diffs are extra Ts in the second, the 'CT' repeat in first, then extra 'GG' in first then extra 'G' in second

        int mismatches = SequenceCompare.compareSequences(
                firstBases.getBytes(), firstBaseQuals, 0, firstBaseQuals.length - 1, firstRepeats,
                secondBases.getBytes(), secondBaseQuals, 0, secondBaseQuals.length - 1, secondRepeats,
                NO_MISMATCH_LIMIT);

        assertEquals(4, mismatches);

        mismatches = SequenceCompare.compareSequences(
                firstBases.getBytes(), firstBaseQuals, 0, firstBaseQuals.length - 1, firstRepeats,
                secondBases.getBytes(), secondBaseQuals, 0, secondBaseQuals.length - 1, secondRepeats,
                NO_MISMATCH_LIMIT);

        assertEquals(4, mismatches);

        // a single SNV at the start of a start of a homopolymer can look like two 1-based INDELs or a repeat count diff
        // but only count this as one mismatch if it can be explained by a single SNV
        firstBases =  "AACGTTTTTAGCTGA";
        secondBases = "AACTTTTTTAGCTGA";

        firstRepeats = findRepeats(firstBases.getBytes());
        secondRepeats = findRepeats(secondBases.getBytes());

        firstBaseQuals = buildDefaultBaseQuals(secondBases.length());
        secondBaseQuals = firstBaseQuals;

        mismatches = SequenceCompare.compareSequences(
                firstBases.getBytes(), firstBaseQuals, 0, firstBaseQuals.length - 1, firstRepeats,
                secondBases.getBytes(), secondBaseQuals, 0, secondBaseQuals.length - 1, secondRepeats,
                NO_MISMATCH_LIMIT);

        assertEquals(1, mismatches);

        mismatches = SequenceCompare.compareSequences(
                firstBases.getBytes(), firstBaseQuals, 0, firstBaseQuals.length - 1, firstRepeats,
                secondBases.getBytes(), secondBaseQuals, 0, secondBaseQuals.length - 1, secondRepeats,
                NO_MISMATCH_LIMIT);

        assertEquals(1, mismatches);
    }

    @Test
    public void testLongerSequenceComparisons()
    {
        //                              10        20            30        40         50            60            70          80        90
        //                    012345678901234567890 1 2345678  901234567890123456789 012  345678 9 012 3 456 7 890123456 7 89012345678901
        String firstBases =  "TTTTTTGTATTAAGTCTAATA C TTTTTTT  AACTTAAGTGTAGATTTTTTT AAA  TGCTCC A TAA C GGT T TTATTTATA C GATTTTTGTCACTG";
        String secondBases = "TTTTTTGTATTAAGTCTAATA G TTTTTTTT AACTTAAGTGTAGATTTTTT  AAAA TGCTCC G TAA T GGT G TTATTTATA T GATTTTTGTCACTG";
        //                    012345678901234567890 1 23456789 01234567890123456789  0123 456789 0 123 4 567 8 901234567 8 90123456789012
        //                              10        20           30        40          50          60             70          80        90

        firstBases = firstBases.replaceAll(" ", "");
        secondBases = secondBases.replaceAll(" ", "");

        List<RepeatInfo> firstRepeats = findRepeats(firstBases.getBytes());
        assertEquals(5, firstRepeats.size());

        byte[] firstBaseQuals = buildDefaultBaseQuals(firstBases.length());

        List<RepeatInfo> secondRepeats = findRepeats(secondBases.getBytes());
        assertEquals(6, secondRepeats.size());

        byte[] secondBaseQuals = buildDefaultBaseQuals(secondBases.length());

        // diffs are extra Ts in the second, the 'CT' repeat in first, then extra 'GG' in first then extra 'G' in second

        int mismatches = SequenceCompare.compareSequences(
                firstBases.getBytes(), firstBaseQuals, 0, firstBaseQuals.length - 1, firstRepeats,
                secondBases.getBytes(), secondBaseQuals, 0, secondBaseQuals.length - 1, secondRepeats,
                NO_MISMATCH_LIMIT);

        assertEquals(8, mismatches);

        mismatches = SequenceCompare.compareSequences(
                firstBases.getBytes(), firstBaseQuals, 0, firstBaseQuals.length - 1, firstRepeats,
                secondBases.getBytes(), secondBaseQuals, 0, secondBaseQuals.length - 1, secondRepeats,
                NO_MISMATCH_LIMIT);

        assertEquals(8, mismatches);
    }

    @Test
    public void testAssemblyRefBasesTrimmed()
    {
        Junction posJunction = new Junction(CHR_1, 60, FORWARD);

        String extensionSequence = "ACGTTCGTAAAAAAGGGGGGACGTACGTCCCC";
        String refBaseSequence = "ACGTAGAGAGAGACGTCCCCACGG";
        String assemblySequence = refBaseSequence + extensionSequence;
        byte[] baseQuals = buildDefaultBaseQuals(assemblySequence.length());

        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 37, assemblySequence, "24M32S");

        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 38, assemblySequence.substring(1), "23M32S");

        JunctionAssembly assembly = new JunctionAssembler(posJunction).processJunction(List.of(read1, read2)).get(0);

        assertEquals("ACGT_AG4_ACGT_C4_ACGG", assembly.refBasesRepeatedTrimmed()); // 4 + 4 + 4 + 2 + 4
        assertEquals(18, assembly.refBaseTrimLength());

        String refBasesTrimmed = buildTrimmedRefBaseSequence(assembly, 12);
        assertEquals("ACGT_AG4_ACGT", refBasesTrimmed);

        Junction negJunction = new Junction(CHR_1, 60, REVERSE);

        assemblySequence = extensionSequence + refBaseSequence;

        assembly = new JunctionAssembly(negJunction, assemblySequence.getBytes(), baseQuals, extensionSequence.length());

        assembly.buildRepeatInfo();

        assertEquals("ACGT_AG4_ACGT_C4_ACGG", assembly.refBasesRepeatedTrimmed());
    }
}

