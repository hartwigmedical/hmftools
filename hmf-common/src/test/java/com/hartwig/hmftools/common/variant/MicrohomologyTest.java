package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MicrohomologyTest {

    @Test
    public void testRepeats() {
        assertHomology("GTAACCAGGAGTGTATT",
                0,
                "CGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAG",
                "CGTAACCAGGAGTGTATT");
    }

    @Test
    public void testSangerExample() {
        assertHomology("TATC", 5, "GCACTGTATCCACTTGATATCATTAT", "GTATCCACTTGA");
        assertHomology("TATC", 6, "GCACTGTATCCACTTGATATCATTAT", "TATCCACTTGAT");
        assertHomology("TATC", 7, "GCACTGTATCCACTTGATATCATTAT", "ATCCACTTGATA");
        assertHomology("TATC", 8, "GCACTGTATCCACTTGATATCATTAT", "TCCACTTGATAT");
        assertHomology("TATC", 9, "GCACTGTATCCACTTGATATCATTAT", "CCACTTGATATC");
    }

    private static void assertHomology(@NotNull final String expected, int position, @NotNull final String sequence,
            @NotNull final String ref) {
        assertEquals(expected, Microhomology.microhomology(position, sequence, ref));
    }
}
