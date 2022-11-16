package com.hartwig.hmftools.cider

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotEquals
import kotlin.test.assertNotNull

class AnchorBlosumSearcherTest
{
    val ighJ1 = VJAnchorTemplate(
        VJGeneType.IGHJ,
        "IGHJ1",
        "01",
        null,
        "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG",
        "TGGGGCCAGGGCACCCTGGTCACCGTCTCC",
        null)

    val ighJ6 = VJAnchorTemplate(
        VJGeneType.IGHJ,
        "IGHJ6","01", null,
        "ATTACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG",
        "TGGGGGCAAGGGACCACGGTCACCGTCTCC",
        null)

    val ighV1_18 = VJAnchorTemplate(
        VJGeneType.IGHV,
        "IGHV1-18",
        "01",
        null,
        "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTTACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACAAACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGA",
        "AGATCTGACGACACGGCCGTGTATTACTGT",
        null)

    val ighV3_7 = VJAnchorTemplate(
        VJGeneType.IGHV,
        "IGHV3-7", "01", null,
        "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGTAGCTATTGGATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTGGCCAACATAAAGCAAGATGGAAGTGAGAAATACTATGTGGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGA",
        "AGAGCCGAGGACACGGCTGTGTATTACTGT",
        null)

    val vjGeneStore = TestCiderGeneDatastore(listOf(ighJ1, ighJ6, ighV1_18, ighV3_7))

    @Test
    fun testFindJAnchor()
    {
        val anchorBlosumSearcher = AnchorBlosumSearcher(vjGeneStore,
            8)

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
        val anchorBlosumMatch: AnchorBlosumMatch? = anchorBlosumSearcher.searchForAnchor(
            testSeq, listOf(VJGeneType.IGHJ), IAnchorBlosumSearcher.Mode.ALLOW_NEG_SIMILARITY)

        assertNotNull(anchorBlosumMatch)
        assertEquals(45, anchorBlosumMatch.anchorStart)
        assertEquals(jAnchorSeq, testSeq.substring(anchorBlosumMatch.anchorStart, anchorBlosumMatch.anchorEnd))
        //assertEquals(15, anchorBlosumMatch.similarityScore)
    }

    @Test
    fun testFindJAnchorPartial()
    {
        val anchorBlosumSearcher = AnchorBlosumSearcher(vjGeneStore, 7)

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
        val anchorBlosumMatch: AnchorBlosumMatch? = anchorBlosumSearcher.searchForAnchor(
            testSeq, listOf(VJGeneType.IGHJ), IAnchorBlosumSearcher.Mode.ALLOW_NEG_SIMILARITY)

        assertNotNull(anchorBlosumMatch)
        assertEquals(45, anchorBlosumMatch.anchorStart)
        assertEquals(jAnchorSeq.dropLast(1), testSeq.substring(anchorBlosumMatch.anchorStart, anchorBlosumMatch.anchorEnd))
        //assertEquals(9, anchorBlosumMatch.similarityScore)
    }

    @Test
    fun testFindVAnchor()
    {
        val anchorBlosumSearcher = AnchorBlosumSearcher(vjGeneStore, 8)

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
        val anchorBlosumMatch: AnchorBlosumMatch? = anchorBlosumSearcher.searchForAnchor(
            testSeq, listOf(VJGeneType.IGHV), IAnchorBlosumSearcher.Mode.ALLOW_NEG_SIMILARITY)

        assertNotNull(anchorBlosumMatch)
        assertEquals(27, anchorBlosumMatch.anchorStart)
        assertEquals(vAnchorSeq, testSeq.substring(anchorBlosumMatch.anchorStart, anchorBlosumMatch.anchorEnd))
        //assertEquals(11, anchorBlosumMatch.similarityScore)
    }

    @Test
    fun testFindVAnchorPartial()
    {
        val anchorBlosumSearcher = AnchorBlosumSearcher(vjGeneStore, 7)

        val vAnchorSeq = "CGAGTCGAAGACACGGCTGTGTATTACTGT".drop(8) // chop off 8
        val jAnchorSeq = "TGGGGCCAAGGGACCACGGTCACCGTCTCC"
        val fullSeq = vAnchorSeq +
                "GCGAGACCGAAATTTTATAGTAATGGCTTGGCGGGTATGGACGTC" +
                jAnchorSeq +
                "CCAGCATAAAACACGATGGAAGTGATA"
        val jAnchorStart = fullSeq.indexOf(jAnchorSeq)

        assertNotEquals(-1, jAnchorStart)

        val testSeq = fullSeq.substring(0, jAnchorStart)

        // we use negative start offset, this means that relative to the seq start, the anchor start
        // could be negative, and would result in a partial match
        //         ----------------------------------- testSeq
        // ++++++++++    template anchor
        // |------|      startOffset
        val anchorBlosumMatch: AnchorBlosumMatch? = anchorBlosumSearcher.searchForAnchor(
            testSeq, listOf(VJGeneType.IGHV), IAnchorBlosumSearcher.Mode.ALLOW_NEG_SIMILARITY)

        assertNotNull(anchorBlosumMatch)
        assertEquals(1, anchorBlosumMatch.anchorStart)
        assertEquals(vAnchorSeq.drop(1), testSeq.substring(anchorBlosumMatch.anchorStart, anchorBlosumMatch.anchorEnd))
        //assertEquals(9, anchorBlosumMatch.similarityScore)
    }
}