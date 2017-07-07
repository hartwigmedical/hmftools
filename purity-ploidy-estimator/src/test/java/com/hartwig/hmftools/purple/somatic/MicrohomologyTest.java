package com.hartwig.hmftools.purple.somatic;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MicrohomologyTest {

    @Test
    public void testSangaExample() {
        assertHomology("TATC", "TATCCACTTGAT", "T", "TATCCACTTGATATCATTAT");
        assertHomology("TATC", "GTATCCACTTGA", "G", "GTATCCACTTGATATCATTAT");
    }

    @Test
    public void testCOLO829() {
        assertHomology("GTAACCAGGAGTGTATT", "CGTAACCAGGAGTGTATT", "C", "CGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGT");
        assertHomology("AA", "AAAG", "A", "AAAGAAAG");
        assertHomology("A", "CA", "C", "CAAAA");
    }

    private void assertHomology(@NotNull final String expected, @NotNull final String ref, @NotNull final String alt, @NotNull final String genome) {
        assertEquals(expected, Microhomology.microhomology(ref, alt, genome));
    }
}
