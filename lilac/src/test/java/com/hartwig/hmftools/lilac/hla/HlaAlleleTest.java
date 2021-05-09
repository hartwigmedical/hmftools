package com.hartwig.hmftools.lilac.hla;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import org.junit.Test;

public class HlaAlleleTest
{
    @Test
    public void testDecode()
    {
        assertContig("A*01:01:01:01", "A*01:01:01:01");
        assertContig("A*01:01:01", "A*01:01:01");
        assertContig("A*26:22", "A*26:22");
        assertContig("A*26", "A*26");
    }

    @Test
    public void testReduce()
    {
        HlaAllele eightDigit = HlaAllele.fromString("A*01:01:01:01");
        assertMatch(HlaAllele.fromString("A*01:01:01"), eightDigit.asSixDigit());
        assertMatch(HlaAllele.fromString("A*01:01"), eightDigit.asFourDigit());
        assertMatch(HlaAllele.fromString("A*01"), eightDigit.asAlleleGroup());
    }

    private void assertContig(String expected, String contig)
    {
        assertEquals(expected, HlaAllele.fromString(contig).toString());
    }

    private void assertMatch(HlaAllele expected, HlaAllele contig)
    {
        assertTrue(expected.toString().equals(contig.toString()));
    }

}
