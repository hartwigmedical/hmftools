package com.hartwig.hmftools.cider

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotEquals
import kotlin.test.assertNotNull

class AnchorBlosumSearcherTest
{
    val ighJ1 = VJGene(
        "IGHJ1*01",
        "IGHJ1",
        "01",
        null,
        "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG",
        "TGGGGCCAGGGCACCCTGGTCACCGTCTCC",
        null)

    val ighJ6 = VJGene(
        "IGHJ6*01", "IGHJ6","01", null,
        "ATTACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG",
        "TGGGGGCAAGGGACCACGGTCACCGTCTCC",
        null)

    val ighV1_18 = VJGene(
        "IGHV1-18*01",
        "IGHV1-18",
        "01",
        null,
        "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA",
        "AGATCTGACGACACGGCCGTGTATTACTGT",
        null)

    val ighV3_7 = VJGene(
        "IGHV3-7*01", "IGHV3-7", "01", null,
        "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGTAGCTATTGGATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTGGCCAACATAAAGCAAGATGGAAGTGAGAAATACTATGTGGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGA",
        "AGAGCCGAGGACACGGCTGTGTATTACTGT",
        null)

    val vjGeneStore = TestVJGeneStore(listOf(ighJ1, ighJ6, ighV1_18, ighV3_7))

    @Test
    fun testFindJAnchor()
    {
        val anchorBlosumSearcher = AnchorBlosumSearcher(vjGeneStore,
            8, true)

        val vAnchorSeq = "CGAGTCGAAGACACGGCTGTGTATTACTGT"
        val jAnchorSeq = "TGGGGCCAAGGGACCACGGTCACCGTCTCC"
        val fullSeq = "TCACTGTCTCTGCAAATGAATGACCTG" +
                vAnchorSeq +
                "GCGAGACCGAAATTTTATAGTAATGGCTTGGCGGGTATGGACGTC" +
                jAnchorSeq +
                "CCAGCATAAAACACGATGGAAGTGATA"
        val vAnchorStart = fullSeq.indexOf(vAnchorSeq)

        assertNotEquals(-1, vAnchorStart)

        val vAnchorEnd = vAnchorStart + vAnchorSeq.length

        val testSeq = fullSeq.substring(vAnchorEnd)
        val anchorBlosumMatch: AnchorBlosumMatch? = anchorBlosumSearcher.searchForAnchor(testSeq, VJGeneType.IGHJ)

        assertNotNull(anchorBlosumMatch)
        assertEquals(45, anchorBlosumMatch.anchorStart)
        assertEquals(jAnchorSeq, testSeq.substring(anchorBlosumMatch.anchorStart, anchorBlosumMatch.anchorEnd))
        assertEquals(15, anchorBlosumMatch.similarityScore)
    }

    @Test
    fun testFindJAnchorPartial()
    {
        val anchorBlosumSearcher = AnchorBlosumSearcher(vjGeneStore, 7, true)

        val vAnchorSeq = "CGAGTCGAAGACACGGCTGTGTATTACTGT"
        val jAnchorSeq = "TGGGGCCAAGGGACCACGGTCACCGTCTCC".dropLast(8) // chop off 8 bases
        val fullSeq = "TCACTGTCTCTGCAAATGAATGACCTG" +
                vAnchorSeq +
                "GCGAGACCGAAATTTTATAGTAATGGCTTGGCGGGTATGGACGTC" +
                jAnchorSeq
        val vAnchorStart = fullSeq.indexOf(vAnchorSeq)

        assertNotEquals(-1, vAnchorStart)

        val vAnchorEnd = vAnchorStart + vAnchorSeq.length

        val testSeq = fullSeq.substring(vAnchorEnd)
        val anchorBlosumMatch: AnchorBlosumMatch? = anchorBlosumSearcher.searchForAnchor(testSeq, VJGeneType.IGHJ)

        assertNotNull(anchorBlosumMatch)
        assertEquals(45, anchorBlosumMatch.anchorStart)
        assertEquals(jAnchorSeq.dropLast(1), testSeq.substring(anchorBlosumMatch.anchorStart, anchorBlosumMatch.anchorEnd))
        assertEquals(9, anchorBlosumMatch.similarityScore)
    }

    @Test
    fun testFindVAnchor()
    {
        val anchorBlosumSearcher = AnchorBlosumSearcher(vjGeneStore, 8, true)

        val vAnchorSeq = "CGAGTCGAAGACACGGCTGTGTATTACTGT"
        val jAnchorSeq = "TGGGGCCAAGGGACCACGGTCACCGTCTCC"
        val fullSeq = "TCACTGTCTCTGCAAATGAATGACCTG" +
                vAnchorSeq +
                "GCGAGACCGAAATTTTATAGTAATGGCTTGGCGGGTATGGACGTC" +
                jAnchorSeq +
                "CCAGCATAAAACACGATGGAAGTGATA"
        val jAnchorStart = fullSeq.indexOf(jAnchorSeq)

        assertNotEquals(-1, jAnchorStart)

        val testSeq = fullSeq.substring(0, jAnchorStart)
        val anchorBlosumMatch: AnchorBlosumMatch? = anchorBlosumSearcher.searchForAnchor(testSeq, VJGeneType.IGHV)

        assertNotNull(anchorBlosumMatch)
        assertEquals(27, anchorBlosumMatch.anchorStart)
        assertEquals(vAnchorSeq, testSeq.substring(anchorBlosumMatch.anchorStart, anchorBlosumMatch.anchorEnd))
        assertEquals(11, anchorBlosumMatch.similarityScore)
    }

    @Test
    fun testFindVAnchorPartial()
    {
        val anchorBlosumSearcher = AnchorBlosumSearcher(vjGeneStore, 7, true)

        val vAnchorSeq = "CGAGTCGAAGACACGGCTGTGTATTACTGT".drop(8) // chop off 8
        val jAnchorSeq = "TGGGGCCAAGGGACCACGGTCACCGTCTCC"
        val fullSeq = vAnchorSeq +
                "GCGAGACCGAAATTTTATAGTAATGGCTTGGCGGGTATGGACGTC" +
                jAnchorSeq +
                "CCAGCATAAAACACGATGGAAGTGATA"
        val jAnchorStart = fullSeq.indexOf(jAnchorSeq)

        assertNotEquals(-1, jAnchorStart)

        val testSeq = fullSeq.substring(0, jAnchorStart)
        val anchorBlosumMatch: AnchorBlosumMatch? = anchorBlosumSearcher.searchForAnchor(testSeq, VJGeneType.IGHV)

        assertNotNull(anchorBlosumMatch)
        assertEquals(1, anchorBlosumMatch.anchorStart)
        assertEquals(vAnchorSeq.drop(1), testSeq.substring(anchorBlosumMatch.anchorStart, anchorBlosumMatch.anchorEnd))
        assertEquals(9, anchorBlosumMatch.similarityScore)
    }
}