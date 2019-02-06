package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class MicrohomologyTest {

    @Test
    public void worksInRepeatSection() {
        String expectedMicrohomology = "GTAACCAGGAGTGTATT";
        String refGenome = "CGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAACCAGGAGTGTATTGTAG";
        String ref = "CGTAACCAGGAGTGTATT";

        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(0, refGenome, ref));
    }

    @Test
    public void worksOnSangerExamples() {
        String expectedMicrohomology = "TATC";
        String refGenome = "GCACTGTATCCACTTGATATCATTAT";

        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(5, refGenome, "GTATCCACTTGA"));
        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(6, refGenome, "TATCCACTTGAT"));
        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(7, refGenome, "ATCCACTTGATA"));
        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(8, refGenome, "TCCACTTGATAT"));
        assertEquals(expectedMicrohomology, Microhomology.microhomologyAtDelete(9, refGenome, "CCACTTGATATC"));
    }

    @Test
    public void testMicrohomologyOnInsert() {
        final String referenceGenomeSequence = "ATGCGATCTTCC";

        assertEquals("TC", Microhomology.microhomologyAtInsert(7,referenceGenomeSequence, "CTC"));
        assertEquals("TT", Microhomology.microhomologyAtInsert(7,referenceGenomeSequence, "CTT"));
        assertEquals("", Microhomology.microhomologyAtInsert(7,referenceGenomeSequence, "CATG"));
        assertEquals("TT", Microhomology.microhomologyAtInsert(7,referenceGenomeSequence, "CTTA"));

        // Note these two examples below are equivalent
        assertEquals("ATC", Microhomology.microhomologyAtInsert(7,referenceGenomeSequence, "CAATC"));
        assertEquals("ATC", Microhomology.microhomologyAtInsert(4,referenceGenomeSequence, "GATCA"));
    }

    @Test
    public void testCommonPrefix() {
        assertEquals("GATC", Microhomology.commonPrefix("GATC", "GATC", 100));
        assertEquals("GATC", Microhomology.commonPrefix("GATCA", "GATC", 100));
        assertEquals("GATC", Microhomology.commonPrefix("GATC", "GATCB", 100));
        assertEquals("GATC", Microhomology.commonPrefix("GATCA", "GATCB", 100));
        assertEquals("GATC", Microhomology.commonPrefix("GATCAAAAAA", "GATCBBB", 100));
    }

    @Test
    public void testCommonSuffix() {
        assertEquals("GATC", Microhomology.commonSuffix("GATC", "GATC"));
        assertEquals("GATC", Microhomology.commonSuffix("AGATC", "GATC"));
        assertEquals("GATC", Microhomology.commonSuffix("GATC", "AGATC"));
        assertEquals("GATC", Microhomology.commonSuffix("TGATC", "AGATC"));
        assertEquals("GATC", Microhomology.commonSuffix("AAAAAAGATC", "BBGATC"));
    }

}
