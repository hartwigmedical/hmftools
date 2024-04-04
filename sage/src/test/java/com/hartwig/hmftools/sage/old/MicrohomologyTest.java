package com.hartwig.hmftools.sage.old;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.sage.old.MicrohomologyContext;
import com.hartwig.hmftools.sage.old.MicrohomologyContextBuilder;

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

        assertEquals(expectedMicrohomology, MicrohomologyContextBuilder.microhomologyAtDelete(0, refGenome, ref));
    }

    @Test
    public void testMicrohomologyOnInsert()
    {
        final String refSequence = "ATGCGATCTTCC";

        assertEquals("TC", MicrohomologyContextBuilder.microhomologyAtInsert(7, refSequence, "CTC"));
        assertEquals("TT", MicrohomologyContextBuilder.microhomologyAtInsert(7, refSequence, "CTT"));
        assertEquals("", MicrohomologyContextBuilder.microhomologyAtInsert(7, refSequence, "CATG"));
        assertEquals("TT", MicrohomologyContextBuilder.microhomologyAtInsert(7, refSequence, "CTTA"));

        // Note these two examples below are equivalent
        assertEquals("ATC", MicrohomologyContextBuilder.microhomologyAtInsert(4, refSequence, "GATCA"));
        assertEquals("ATC", MicrohomologyContextBuilder.microhomologyAtInsert(5, refSequence, "ATCAA"));
        assertEquals("ATC", MicrohomologyContextBuilder.microhomologyAtInsert(6, refSequence, "TCAAT"));
        assertEquals("ATC", MicrohomologyContextBuilder.microhomologyAtInsert(7, refSequence, "CAATC"));
    }

    @Test
    public void testMicrohomologyOnInsertWithReadSequence()
    {
        final String readSequence = "ATGCGATCAATCTTCC";
        assertEquals("ATC", MicrohomologyContextBuilder.microhomologyAtInsert(4, 5, readSequence.getBytes()).toString());
        assertEquals("ATC", MicrohomologyContextBuilder.microhomologyAtInsert(5, 5, readSequence.getBytes()).toString());
        assertEquals("ATC", MicrohomologyContextBuilder.microhomologyAtInsert(6, 5, readSequence.getBytes()).toString());
        assertEquals("ATC", MicrohomologyContextBuilder.microhomologyAtInsert(7, 5, readSequence.getBytes()).toString());
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
        assertEquals("GATCAA", new String(MicrohomologyContextBuilder.reconstructDeletedSequence(0, "GTCAA".getBytes(), "GA")));
        assertEquals("GATCAA", new String(MicrohomologyContextBuilder.reconstructDeletedSequence(1, "GACAA".getBytes(), "AT")));
        assertEquals("GATCAA", new String(MicrohomologyContextBuilder.reconstructDeletedSequence(2, "GATAA".getBytes(), "TC")));
        assertEquals("GATCAA", new String(MicrohomologyContextBuilder.reconstructDeletedSequence(3, "GATCA".getBytes(), "CA")));
        assertEquals("GATCAA", new String(MicrohomologyContextBuilder.reconstructDeletedSequence(4, "GATCA".getBytes(), "AA")));
    }

    @Test
    public void testReconstructDeletedSequenceWithDelCombinedWithSnv()
    {
        assertEquals("GATCAA", new String(MicrohomologyContextBuilder.reconstructDeletedSequence(1, "GACAA".getBytes(), "TT")));
        assertEquals("GATCAA", new String(MicrohomologyContextBuilder.reconstructDeletedSequence(1, "GACAA".getBytes(), "AT")));
    }

    private static void assertDeleteExpandRepeats(@NotNull String expected, int position, int altLength, @NotNull String readSequence)
    {
        MicrohomologyContext context = MicrohomologyContextBuilder.microhomologyAtDelete(position, altLength, readSequence.getBytes());
        MicrohomologyContext context2 = MicrohomologyContextBuilder.expandMicrohomologyRepeats(context);
        assertEquals(expected, context2.toString());
    }

    private static void assertInsertExpandRepeats(@NotNull String expected, int position, int altLength, @NotNull String readSequence)
    {
        MicrohomologyContext context = MicrohomologyContextBuilder.microhomologyAtInsert(position, altLength, readSequence.getBytes());
        MicrohomologyContext context2 = MicrohomologyContextBuilder.expandMicrohomologyRepeats(context);
        assertEquals(expected, context2.toString());
    }

    private static void assertInsert(@NotNull String expected, int position, int altLength, @NotNull String readSequence)
    {
        MicrohomologyContext context = MicrohomologyContextBuilder.microhomologyAtInsert(position, altLength, readSequence.getBytes());
        assertEquals(expected, context.toString());
    }
}
