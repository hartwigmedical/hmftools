package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
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

import java.util.List;

import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.assembly.read.Read;

import org.junit.Test;

public class SequenceTest
{
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
                NO_MISMATCH_LIMIT, false, true);

        assertEquals(2, mismatches);
    }

        @Test
    public void testSequenceComparisons()
    {
        //                                   10          20            30
        //                    012345     6789012345 678 90 12345678   9012
        String firstBases =  "ATTTTT     AACTCTCTCT AAA GG CTGACGTA   TTCC";

        String secondBases = "ATTTTTTTTT AACTCTCT   AAA    CTGACGTA G TTCC";
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
                NO_MISMATCH_LIMIT, false, true);

        assertEquals(4, mismatches);

        // a single SNV at the start of a start of a homopolymer can look like two 1-based INDELs or a repeat count diff
        // but only count this as one mismatch if it can be explained by a single SNV
        firstBases =  "AACGTTTTTAGCTGA";
        secondBases = "AACTTTTTTAGCTGA";

        firstRepeats = findRepeats(secondBases.getBytes());
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
                NO_MISMATCH_LIMIT, false, true);

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
                NO_MISMATCH_LIMIT, false, true);

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
        Read read1b = cloneRead(read1, READ_ID_GENERATOR.nextId());

        JunctionAssembly assembly = new JunctionAssembler(posJunction).processJunction(List.of(read1, read1b)).get(0);

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

