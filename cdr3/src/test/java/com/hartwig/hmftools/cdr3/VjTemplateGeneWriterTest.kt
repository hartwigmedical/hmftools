package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.common.genome.region.Strand
import kotlin.test.*

class VjTemplateGeneWriterTest
{
    @Test
    fun testCalcPositionMismatch()
    {
        val seq = "ABCDE"
        val refSeq = "CDEFG"

        val (startShift, numMisMatch) = VjTemplateGeneWriter.calcPositionMismatch(seq, refSeq, 3, 3)

        assertEquals(-2, startShift)
        assertEquals(0, numMisMatch)
    }

    /*
    @Test
    fun testCalcPositionMismatchLong()
    {
        val seq = "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTTCAGCCTCTGGATTCACCTTCAGTAGCTATGCTATGCACTGGGTCCGCCAGGCTC" +
                "CAGGGAAGGGACTGGAATATGTTTCAGCTATTAGTAGTAATGGGGGTAGCACATACTACGCAGACTCAGTGAAGGGCAGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATGTCC" +
                "AAATGAGCAGTCTGAGAGCTGAGGACACGGCTGTGTATTACTGTGTGAAAGA"
        val refSeq = "AGCTCTGGGAGAGGAGCCCCCGCCCTGGGATTCCCAGGTGTTTTCATTTGGTGATCAGCACTGAACACAGAAGAGTCATGATGGAGTTTGGGCTGAGCTGGGTTTTCCTTGTTGCTAT" +
                "TTTTAAAGGTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGAAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATGCTATGCA" +
                "CTGGGTCCGCCAGGCTCCAGGGAAGGGACTGGAATATGTTTCAGCTATTAGTAGTAATGGGGGTAGCACATATTATGCAGACTCTGTGAAGGGCAGATTCACCATCTCCAGAGACAATTCCAA" +
                "GAACACGCTGTATCTTCAAATGGGCAGCCTGAGAGCTGAGGACATGGCTGTGTATTACTGTGCGAGA"

        val (startShift, numMisMatch) = ImgtAnchorFinder.calcPositionMismatch(seq, refSeq)

        assertEquals(-2, startShift)
        assertEquals(0, numMisMatch)
    }*/

    @Test
    fun testCorrectGeneLocationShort()
    {
        val seq = "ABCDE"
        val refSeq = "CDEFGH"

        var refGeneLoc = GeneLocation("1", 100, 105, Strand.FORWARD)
        var correctedGeneLoc = VjTemplateGeneWriter.correctGeneLocation(seq, refSeq, refGeneLoc, 3, 3)

        // 2 more bases at the start
        assertEquals(98, correctedGeneLoc.posStart)
        // 3 less bases at the end
        assertEquals(102, correctedGeneLoc.posEnd)

        // now test with reverse strand
        refGeneLoc = GeneLocation("1", 100, 105, Strand.REVERSE)
        correctedGeneLoc = VjTemplateGeneWriter.correctGeneLocation(seq, refSeq, refGeneLoc, 3, 3)

        // 3 less bases at the start
        assertEquals(103, correctedGeneLoc.posStart)
        // 2 more bases at the end
        assertEquals(107, correctedGeneLoc.posEnd)
    }

    @Test
    fun testCorrectGeneLocationSimple()
    {
        val seq = "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG"
        val refSeq = "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG"

        var refGeneLoc = GeneLocation("1", 100, 151, Strand.FORWARD)
        var correctedGeneLoc = VjTemplateGeneWriter.correctGeneLocation(seq, refSeq, refGeneLoc, 20, 30)

        assertEquals(refGeneLoc.posStart, correctedGeneLoc.posStart)
        // 3 less bases at the end
        assertEquals(refGeneLoc.posEnd, correctedGeneLoc.posEnd)
    }

    @Test
    fun testCorrectGeneLocation1()
    {
        val seq = "CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCATCTGCGCTGTCTCTGGTGACTCCATCAGCAGTGGTAACTGGTGAATCTGGGTCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCCATCATAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAATCACCATGTCCGTAGACACGTCCAAGAACCAGTTCTACCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCCGTGTATTACTGTGCGAGATA"
        val refSeq = "ATGAAACACCTGTGGTTCTTCCTCCTGCTGGTGGCAGCTCCCAGATGGGTCCTGTCCCAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCATCTGCGCTGTCTCTGGTGACTCCATCAGCAGTGGTAACTGGTGAATCTGGGTCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCCATCATAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAATCACCATGTCCGTAGACACGTCCAAGAACCAGTTCTACCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCCGTGTATTACTGTGCGAGATACACAGTGAGGGGAGG"

        var refGeneLoc = GeneLocation("1", 1, 300, Strand.FORWARD)
        var correctedGeneLoc = VjTemplateGeneWriter.correctGeneLocation(seq, refSeq, refGeneLoc, 20, 30)

        assertEquals(58, correctedGeneLoc.posStart)
        // 3 less bases at the end
        assertEquals(285, correctedGeneLoc.posEnd)
    }
}