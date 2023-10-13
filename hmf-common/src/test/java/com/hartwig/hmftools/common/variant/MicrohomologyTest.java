package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MicrohomologyTest
{
    @Test
    public void worksInRepeatSection()
    {
        String expectedMicrohomology = "GTAACCAGGAGTGTATT";
        String refGenome = "CGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAG";
        String ref = "CGTAACCAGGAGTGTATT";

        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(0, refGenome, ref));
    }

    @Test
    public void worksOnSangerExamples()
    {
        String expectedMicrohomology = "TATC";
        String refGenome = "GCACTGTATCCACTTGATATCATTAT";

        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(5, refGenome, "GTATCCACTTGA"));
        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(6, refGenome, "TATCCACTTGAT"));
        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(7, refGenome, "ATCCACTTGATA"));
        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(8, refGenome, "TCCACTTGATAT"));
        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(9, refGenome, "CCACTTGATATC"));
    }

    @Test
    public void testMicrohomologyOnInsert()
    {
        final String refSequence = "ATGCGATCTTCC";

        assertEquals("TC", Microhomology.microhomologyAtInsert(7, refSequence, "CTC"));
        assertEquals("TT", Microhomology.microhomologyAtInsert(7, refSequence, "CTT"));
        assertEquals("", Microhomology.microhomologyAtInsert(7, refSequence, "CATG"));
        assertEquals("TT", Microhomology.microhomologyAtInsert(7, refSequence, "CTTA"));

        // Note these two examples below are equivalent
        assertEquals("ATC", Microhomology.microhomologyAtInsert(4, refSequence, "GATCA"));
        assertEquals("ATC", Microhomology.microhomologyAtInsert(5, refSequence, "ATCAA"));
        assertEquals("ATC", Microhomology.microhomologyAtInsert(6, refSequence, "TCAAT"));
        assertEquals("ATC", Microhomology.microhomologyAtInsert(7, refSequence, "CAATC"));
    }

    @Test
    public void testMicrohomologyOnInsertWithReadSequence()
    {
        final String readSequence = "ATGCGATCAATCTTCC";
        assertEquals("ATC", Microhomology.microhomologyAtInsert(4, 5, readSequence.getBytes()).toString());
        assertEquals("ATC", Microhomology.microhomologyAtInsert(5, 5, readSequence.getBytes()).toString());
        assertEquals("ATC", Microhomology.microhomologyAtInsert(6, 5, readSequence.getBytes()).toString());
        assertEquals("ATC", Microhomology.microhomologyAtInsert(7, 5, readSequence.getBytes()).toString());
    }

    @Test
    public void testMicrohomologyInsertInRepeat()
    {
        final String readSequence = "ATTTTGTTTGTTTGA";

        assertInsert("", 0, 5, readSequence);
        assertInsert("TTTG", 1, 5, readSequence);
        assertInsert("TTTG", 2, 5, readSequence);
        assertInsert("TTTG", 3, 5, readSequence);
        assertInsert("TTTG", 4, 5, readSequence);
        assertInsert("TTTG", 5, 5, readSequence);
        assertInsert("TTTG", 6, 5, readSequence);
        assertInsert("TTTG", 7, 5, readSequence);
        assertInsert("TTTG", 8, 5, readSequence);
        assertInsert("TTTG", 9, 5, readSequence);
        assertInsert("", 10, 5, readSequence);
    }

    @Test
    public void testMicrohomologyInsertWithExpandedRepeats()
    {
        final String readSequence = "ATTTTGTTTGTTTGA";

        assertInsertExpandRepeats("", 0, 5, readSequence);
        assertInsertExpandRepeats("TTTGTTTGTTTG", 1, 5, readSequence);
        assertInsertExpandRepeats("TTTGTTTGTTTG", 2, 5, readSequence);
        assertInsertExpandRepeats("TTTGTTTGTTTG", 3, 5, readSequence);
        assertInsertExpandRepeats("TTTGTTTGTTTG", 4, 5, readSequence);
        assertInsertExpandRepeats("TTTGTTTGTTTG", 5, 5, readSequence);
        assertInsertExpandRepeats("TTTGTTTGTTTG", 6, 5, readSequence);
        assertInsertExpandRepeats("TTTGTTTGTTTG", 7, 5, readSequence);
        assertInsertExpandRepeats("TTTGTTTGTTTG", 8, 5, readSequence);
        assertInsertExpandRepeats("TTTGTTTGTTTG", 9, 5, readSequence);
        assertInsertExpandRepeats("", 10, 5, readSequence);
    }

    @Test
    public void testMicrohomologyDeleteWithExpandedRepeats()
    {
        final String refSequence = "ATTTTGTTTGTTTGA";

        assertDeleteExpandRepeats("", 0, 5, refSequence);
        assertDeleteExpandRepeats("TTTGTTTGTTTG", 1, 5, refSequence);
        assertDeleteExpandRepeats("TTTGTTTGTTTG", 2, 5, refSequence);
        assertDeleteExpandRepeats("TTTGTTTGTTTG", 3, 5, refSequence);
        assertDeleteExpandRepeats("TTTGTTTGTTTG", 4, 5, refSequence);
        assertDeleteExpandRepeats("TTTGTTTGTTTG", 5, 5, refSequence);
        assertDeleteExpandRepeats("TTTGTTTGTTTG", 6, 5, refSequence);
        assertDeleteExpandRepeats("TTTGTTTGTTTG", 7, 5, refSequence);
        assertDeleteExpandRepeats("TTTGTTTGTTTG", 8, 5, refSequence);
        assertDeleteExpandRepeats("TTTGTTTGTTTG", 9, 5, refSequence);
        assertDeleteExpandRepeats("", 10, 5, refSequence);
    }

    @Test
    public void testReconstructDeletedSequence()
    {
        assertEquals("GATCAA", new String(Microhomology.reconstructDeletedSequence(0, "GTCAA".getBytes(), "GA")));
        assertEquals("GATCAA", new String(Microhomology.reconstructDeletedSequence(1, "GACAA".getBytes(), "AT")));
        assertEquals("GATCAA", new String(Microhomology.reconstructDeletedSequence(2, "GATAA".getBytes(), "TC")));
        assertEquals("GATCAA", new String(Microhomology.reconstructDeletedSequence(3, "GATCA".getBytes(), "CA")));
        assertEquals("GATCAA", new String(Microhomology.reconstructDeletedSequence(4, "GATCA".getBytes(), "AA")));
    }

    @Test
    public void testReconstructDeletedSequenceWithDelCombinedWithSnv()
    {
        assertEquals("GATCAA", new String(Microhomology.reconstructDeletedSequence(1, "GACAA".getBytes(), "TT")));
        assertEquals("GATCAA", new String(Microhomology.reconstructDeletedSequence(1, "GACAA".getBytes(), "AT")));
    }

    private static void assertDeleteExpandRepeats(@NotNull String expected, int position, int altLength, @NotNull String readSequence)
    {
        MicrohomologyContext context = Microhomology.microhomologyAtDelete(position, altLength, readSequence.getBytes());
        MicrohomologyContext context2 = Microhomology.expandMicrohomologyRepeats(context);
        assertEquals(expected, context2.toString());
    }

    private static void assertInsertExpandRepeats(@NotNull String expected, int position, int altLength, @NotNull String readSequence)
    {
        MicrohomologyContext context = Microhomology.microhomologyAtInsert(position, altLength, readSequence.getBytes());
        MicrohomologyContext context2 = Microhomology.expandMicrohomologyRepeats(context);
        assertEquals(expected, context2.toString());
    }

    private static void assertInsert(@NotNull String expected, int position, int altLength, @NotNull String readSequence)
    {
        MicrohomologyContext context = Microhomology.microhomologyAtInsert(position, altLength, readSequence.getBytes());
        assertEquals(expected, context.toString());
    }
}
