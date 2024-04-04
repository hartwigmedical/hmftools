package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.RepeatInfo.findMaxRepeat;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class VariantReadContextTest
{
    @Test
    public void testRepeats()
    {
        //              01234567890123456789
        String bases = "ACGTACGTTTTTTGTACGT";

        int maxLength = 5;
        int minCount = 3;

        int indexEnd = bases.length() - 1;

        RepeatInfo repeat = findMaxRepeat(bases.getBytes(), 0, indexEnd, maxLength, minCount, false, -1);

        assertNotNull(repeat);
        assertEquals("T", repeat.Bases);
        assertEquals(7, repeat.Index);
        assertEquals(6, repeat.Count);

        repeat = findMaxRepeat(bases.getBytes(), 10, indexEnd, maxLength, minCount, true, -1);

        assertNotNull(repeat);
        assertEquals("T", repeat.Bases);
        assertEquals(7, repeat.Index);
        assertEquals(6, repeat.Count);

        // longer repeats

        //       01234567890123456789
        bases = "ACGTACGAGAGAGGTACGT";

        repeat = findMaxRepeat(bases.getBytes(), 0, indexEnd, maxLength, minCount, true, -1);

        assertNotNull(repeat);
        assertEquals("GA", repeat.Bases);
        assertEquals(6, repeat.Index);
        assertEquals(3, repeat.Count);

        //       0123456789012345678901234567890123456789
        bases = "ACGTACACGGAAAGGAAAGGAAAGGAAAGGAAAGTACGT";

        repeat = findMaxRepeat(bases.getBytes(), 0, indexEnd, maxLength, minCount, true, -1);

        assertNotNull(repeat);
        assertEquals("GGAAA", repeat.Bases);
        assertEquals(8, repeat.Index);
        assertEquals(5, repeat.Count);

        // finding and extending the start
        repeat = findMaxRepeat(bases.getBytes(), 17, indexEnd, maxLength, minCount, true, -1);

        assertNotNull(repeat);
        assertEquals("GGAAA", repeat.Bases);
        assertEquals(8, repeat.Index);
        assertEquals(5, repeat.Count);

        // longest when there are lots to choose from

        //                 10        20        30        40        50
        //       012345678901234567890123456789012345678901234567890123456789
        bases = "AAAACGCGCGCGTTTTTGCTTTGCTTTGCTTTAACGTACGTACGTACGTCCDDCCCDDD";

        repeat = findMaxRepeat(bases.getBytes(), 0, bases.length() - 1, maxLength, minCount, false, 20);

        assertNotNull(repeat);
        assertEquals("TTTGC", repeat.Bases);
        assertEquals(14, repeat.Index);
        assertEquals(3, repeat.Count);

        // checks that the repeat crosses the required index
        repeat = findMaxRepeat(bases.getBytes(), 0, bases.length() - 1, maxLength, minCount, false, 40);

        assertNotNull(repeat);
        assertEquals("ACGT", repeat.Bases);
        assertEquals(33, repeat.Index);
        assertEquals(4, repeat.Count);

        repeat = findMaxRepeat(bases.getBytes(), 0, bases.length() - 1, maxLength, minCount, false, 8);

        assertNotNull(repeat);
        assertEquals("CG", repeat.Bases);
        assertEquals(4, repeat.Index);
        assertEquals(4, repeat.Count);
    }

    private static final String REF_BASES =
        //             10        20        30        40        50        60        70        80        90
        //   0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
            "CGCAATATTCGGGTGGGAGTGACCCGATTTTCCAGGTGCGTTCGTCACCGCTGTCTGTGACTCGGAAAGGGAACTCCCTGACCCCTTGCGCTTCCCAGGT"
          + "GAGGCAATGCCTCGCCCTGCTTCGGCTCGCGCACAGTGCGCGCACACACTGGCCTGCGCCCACTGTCTGGCACTCCCTAGTGAGATGAACCCGGTACCTC";

    private static final RefSequence REF_SEQUENCE = new RefSequence(0, REF_BASES.getBytes()); // note zero-based to line up with indices

    @Test
    public void testSnvMnvVariantReadContext()
    {
        SimpleVariant var = createSimpleVariant(50, "C", "A");

        String readBases = REF_BASES.substring(30, 50) + "A" + REF_BASES.substring(51, 70);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "40M";
        SAMRecord read = buildSamRecord(30, readCigar, readBases, baseQuals);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(5);

        VariantReadContext readContext = builder.createMnvContext(var, read, 20, REF_SEQUENCE);

        assertEquals(43, readContext.AlignmentStart);
        assertEquals(57, readContext.AlignmentEnd);
        assertEquals(5, readContext.CoreIndexStart);
        assertEquals(7, readContext.VarReadIndex);
        assertEquals(9, readContext.CoreIndexEnd);
        assertEquals("GTCACCGCTGTCTGT", readContext.refBases());
        assertEquals("GTCACCGATGTCTGT", readContext.readBases());
        assertEquals("GCT", readContext.trinucleotideStr());
        assertEquals(5, readContext.coreLength());
        assertEquals(5, readContext.leftFlankLength());
        assertEquals(5, readContext.rightFlankLength());
        assertEquals("CGATG", readContext.coreStr());
        assertEquals("GTCAC", readContext.leftFlankStr());
        assertEquals("TCTGT", readContext.rightFlankStr());
    }
}
