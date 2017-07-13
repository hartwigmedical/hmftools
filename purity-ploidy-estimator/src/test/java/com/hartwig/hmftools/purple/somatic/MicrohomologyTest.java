package com.hartwig.hmftools.purple.somatic;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MicrohomologyTest {

    @Test
    public void testRepeats() {
        assertHomology("GTAACCAGGAGTGTATT",
                0,
                "CGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAG",
                "CGTAACCAGGAGTGTATT",
                "C");
    }

    @Test
    public void testSangerExample() {
        assertHomology("TATC", 5, "GCACTGTATCCACTTGATATCATTAT", "GTATCCACTTGA", "G");
        assertHomology("TATC", 6, "GCACTGTATCCACTTGATATCATTAT", "TATCCACTTGAT", "T");
        assertHomology("TATC", 7, "GCACTGTATCCACTTGATATCATTAT", "ATCCACTTGATA", "A");
        assertHomology("TATC", 8, "GCACTGTATCCACTTGATATCATTAT", "TCCACTTGATAT", "T");
        assertHomology("TATC", 9, "GCACTGTATCCACTTGATATCATTAT", "CCACTTGATATC", "C");
    }

    private void assertHomology(@NotNull final String expected, int position, @NotNull final String sequence, @NotNull final String ref,
            @NotNull final String alt) {
        assertEquals(expected, Microhomology.microhomology(position, sequence, ref, alt));
    }
}
