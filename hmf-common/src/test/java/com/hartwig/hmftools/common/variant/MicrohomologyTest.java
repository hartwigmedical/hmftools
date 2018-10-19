package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MicrohomologyTest {

    @Test
    public void worksInRepeatSection() {
        String expectedMicrohomology = "GTAACCAGGAGTGTATT";
        String refGenome = "CGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAG";
        String ref = "CGTAACCAGGAGTGTATT";

        assertEquals(expectedMicrohomology, Microhomology.microhomology(0, refGenome, ref));
    }

    @Test
    public void worksOnSangerExamples() {
        String expectedMicrohomology = "TATC";
        String refGenome = "GCACTGTATCCACTTGATATCATTAT";

        assertEquals(expectedMicrohomology, Microhomology.microhomology(5, refGenome, "GTATCCACTTGA"));
        assertEquals(expectedMicrohomology, Microhomology.microhomology(6, refGenome, "TATCCACTTGAT"));
        assertEquals(expectedMicrohomology, Microhomology.microhomology(7, refGenome, "ATCCACTTGATA"));
        assertEquals(expectedMicrohomology, Microhomology.microhomology(8, refGenome, "TCCACTTGATAT"));
        assertEquals(expectedMicrohomology, Microhomology.microhomology(9, refGenome, "CCACTTGATATC"));
    }
}
