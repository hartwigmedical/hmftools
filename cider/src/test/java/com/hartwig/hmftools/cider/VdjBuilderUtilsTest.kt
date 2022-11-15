package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.cider.layout.ReadLayoutBuilderTest
import com.hartwig.hmftools.cider.layout.TestLayoutRead
import htsjdk.samtools.SAMUtils
import org.junit.Before
import kotlin.test.*

class VdjBuilderUtilsTest
{
    @Before
    fun setUp()
    {
        //org.apache.logging.log4j.core.config.Configurator.setRootLevel(org.apache.logging.log4j.Level.TRACE)
    }

    // test the function that identifies overlaps
    @Test
    fun testFindVjOverlap1()
    {
        val seq1 = "GACTCAGTCCAGGAAAGTATT"
        val seq2 =       "GTCCAGGAAAGTATTATACTCA"

        var overlap = VdjBuilderUtils.findSequenceOverlap(seq1, seq2, 10)
        assertNotNull(overlap)
        assertEquals(0, overlap.seq1Offset)
        assertEquals(6, overlap.seq2Offset)
        assertEquals(15, overlap.overlapBases)
        assertEquals(15, overlap.highQualMatchBases)

        // try the other way round
        overlap = VdjBuilderUtils.findSequenceOverlap(seq2, seq1, 10)
        assertNotNull(overlap)
        assertEquals(6, overlap.seq1Offset)
        assertEquals(0, overlap.seq2Offset)
        assertEquals(15, overlap.overlapBases)
        assertEquals(15, overlap.highQualMatchBases)

        // set higher threshold
        overlap = VdjBuilderUtils.findSequenceOverlap(seq1, seq2, 16)
        assertNull(overlap)

        overlap = VdjBuilderUtils.findSequenceOverlap(seq2, seq1, 16)
        assertNull(overlap)
    }

    // test one long one short
    @Test
    fun testFindVjOverlapLongShort()
    {
        val seq1 = "GACTCAGTCCAGGAAAGTATTATACTCA"
        val seq2 =       "GTCCAGGAAAGTATT"

        var overlap = VdjBuilderUtils.findSequenceOverlap(seq1, seq2, 10)
        assertNotNull(overlap)
        assertEquals(0, overlap.seq1Offset)
        assertEquals(6, overlap.seq2Offset)
        assertEquals(15, overlap.overlapBases)

        // try the other way round
        overlap = VdjBuilderUtils.findSequenceOverlap(seq2, seq1, 10)
        assertNotNull(overlap)
        assertEquals(6, overlap.seq1Offset)
        assertEquals(0, overlap.seq2Offset)
        assertEquals(15, overlap.overlapBases)
    }

    @Test
    fun testGetCompareRanges()
    {
        // vdj1: GATG-CGA-ATACC
        // vdj2:   TG-CGA-ATACCAG

        val vdj1: VDJSequence = TestUtils.createVDJ("GATGCGAATACC", 4, 7)
        val vdj2: VDJSequence = TestUtils.createVDJ("TGCGAATACCAG", 2, 5)

        val (vdj1Range: IntRange, vdj2Range: IntRange) = VdjBuilderUtils.getCompareRanges(vdj1, vdj2)

        assertEquals(2, vdj1Range.first)
        assertEquals(11, vdj1Range.last)

        assertEquals(0, vdj2Range.first)
        assertEquals(9, vdj2Range.last)
    }

    @Test
    fun testGetCompareRangesVOnly()
    {
        // vdj1: GATG-CGAATACC
        // vdj2:   TG-CGAATACCAG

        val vdj1: VDJSequence = TestUtils.createVDJ("GATGCGAATACC", 4, null)
        val vdj2: VDJSequence = TestUtils.createVDJ("TGCGAATACCAG", 2, null)

        val (vdj1Range: IntRange, vdj2Range: IntRange) = VdjBuilderUtils.getCompareRanges(vdj1, vdj2)

        assertEquals(2, vdj1Range.first)
        assertEquals(11, vdj1Range.last)

        assertEquals(0, vdj2Range.first)
        assertEquals(9, vdj2Range.last)
    }

    // for one sided sequences, the side without anchor is capped to 45 bases
    @Test
    fun testGetCompareRangesVOnlyLong()
    {
        val vdj1: VDJSequence = TestUtils.createVDJ("T".repeat(100), 4, null)
        val vdj2: VDJSequence = TestUtils.createVDJ("T".repeat(100), 2, null)

        val (vdj1Range: IntRange, vdj2Range: IntRange) = VdjBuilderUtils.getCompareRanges(vdj1, vdj2)

        assertEquals(2, vdj1Range.first)
        assertEquals(4 + CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES - 1, vdj1Range.last)

        assertEquals(0, vdj2Range.first)
        assertEquals(2 + CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES - 1, vdj2Range.last)
    }

    @Test
    fun testGetCompareRangesJOnly()
    {
        // vdj1:     TGCGAAT-ACCAG
        // vdj2:   GATGCGAAT-ACC

        val vdj1: VDJSequence = TestUtils.createVDJ("TGCGAATACCAG", 7, null)
        val vdj2: VDJSequence = TestUtils.createVDJ("GATGCGAATACC", 9, null)

        val (vdj1Range: IntRange, vdj2Range: IntRange) = VdjBuilderUtils.getCompareRanges(vdj1, vdj2)

        assertEquals(0, vdj1Range.first)
        assertEquals(9, vdj1Range.last)

        assertEquals(2, vdj2Range.first)
        assertEquals(11, vdj2Range.last)
    }

    // for one sided sequences, the side without anchor is capped to 60 bases
    @Test
    fun testGetCompareRangesJOnlyLong()
    {
        val vdj1: VDJSequence = TestUtils.createVDJ("T".repeat(100), null, 90)
        val vdj2: VDJSequence = TestUtils.createVDJ("T".repeat(100), null, 95)

        val (vdj1Range: IntRange, vdj2Range: IntRange) = VdjBuilderUtils.getCompareRanges(vdj1, vdj2)

        assertEquals(90 - CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES, vdj1Range.first)
        assertEquals(94, vdj1Range.last)

        assertEquals(95 - CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES, vdj2Range.first)
        assertEquals(99, vdj2Range.last)
    }

    @Test
    fun testMergeVDJsSimple()
    {
        // create two VDJs, both simple read match, exactly the same
        val seq = "ATGCTGGTGT"

        val layout1 = ReadLayout()

        // we are aligned at the T
        val baseQual1 = SAMUtils.fastqToPhred("FF:FFFF:FF") // F is 37, : is 25
        val read1 = TestLayoutRead("read1", ReadKey("read1", true), seq, baseQual1, 4)
        layout1.addRead(read1, ReadLayoutBuilderTest.MIN_BASE_QUALITY)

        val layout2 = ReadLayout()
        val baseQual2 = SAMUtils.fastqToPhred("FFFF::FFFF") // F is 37, : is 25
        val read2 = TestLayoutRead("read2", ReadKey("read2", true), seq, baseQual2, 4)
        layout2.addRead(read2, ReadLayoutBuilderTest.MIN_BASE_QUALITY)

        val vdj1 = TestUtils.createVDJ(layout1, 3, 7, 0, 10)
        val vdj2 = TestUtils.createVDJ(layout2, 3, 7, 0, 10)

        val vdjCombine = VdjBuilderUtils.mergeVDJs(vdj1, vdj2, VDJSequenceBuilderTest.MIN_BASE_QUALITY)

        assertEquals(seq, vdjCombine.sequence)
        assertEquals("2212112122", vdjCombine.supportString) // due to low base qual in some reads
    }

    @Test
    fun testMergeVDJsLongShort()
    {
        // layout1:        CCTCA-GGTG-AT
        // VDJ1:             TCA-GGTG-AT
        // layout2:            A-GGTG-ATAA
        // VDJ2:               A-GGTG-ATA
        // merge layout:   CCTCA-GGTG-ATAA
        // merge VDJ:        TCA-GGTG-ATA

        val layout1 = TestUtils.createLayout("CCTCA-GGTG-AT".replace("-", ""), 5)
        val vdj1 = TestUtils.createVDJ(layout1, 3, 7, 2, 11)
        assertEquals("TCA-GGTG-AT", vdj1.sequenceFormatted)

        val layout2 = TestUtils.createLayout("A-GGTG-ATAA".replace("-", ""), 3)
        val vdj2 = TestUtils.createVDJ(layout2, 1, 5, 0, 8)
        assertEquals("A-GGTG-ATA", vdj2.sequenceFormatted)

        var vdjCombine = VdjBuilderUtils.mergeVDJs(vdj1, vdj2, VDJSequenceBuilderTest.MIN_BASE_QUALITY)

        assertNotNull(vdjCombine.vAnchor)
        assertNotNull(vdjCombine.jAnchor)
        assertEquals("CCTCA-GGTG-ATAA".replace("-", ""), vdjCombine.layout.consensusSequence())
        assertEquals("TCA-GGTG-ATA", vdjCombine.sequenceFormatted)
        assertEquals(3, vdjCombine.vAnchor!!.anchorBoundary)
        assertEquals(7, vdjCombine.jAnchor!!.anchorBoundary)

        // try to merge it the other way, result should be the same
        vdjCombine = VdjBuilderUtils.mergeVDJs(vdj2, vdj1, VDJSequenceBuilderTest.MIN_BASE_QUALITY)

        assertNotNull(vdjCombine.vAnchor)
        assertNotNull(vdjCombine.jAnchor)
        assertEquals("CCTCA-GGTG-ATAA".replace("-", ""), vdjCombine.layout.consensusSequence())
        assertEquals("TCA-GGTG-ATA", vdjCombine.sequenceFormatted)
        assertEquals(3, vdjCombine.vAnchor!!.anchorBoundary)
        assertEquals(7, vdjCombine.jAnchor!!.anchorBoundary)
    }

    @Test
    fun testMergeVDJsVOnly()
    {
        // merge two V only VDJs
        // layout1:        CCTCA-GGTG
        // VDJ1:             TCA-GGTG
        // layout2:            A-GGTGATAA
        // VDJ2:               A-GGTGAT
        // merge layout:   CCTCA-GGTGATAA
        // merge VDJ:        TCA-GGTGAT

        val layout1 = TestUtils.createLayout("CCTCA-GGTG".replace("-", ""), 5)
        val vdj1 = TestUtils.createVDJ(layout1, 3, null, 2, 9)
        assertEquals("TCA-GGTG", vdj1.sequenceFormatted)

        val layout2 = TestUtils.createLayout("A-GGTGATAA".replace("-", ""), 3)
        val vdj2 = TestUtils.createVDJ(layout2, 1, null, 0, 7)
        assertEquals("A-GGTGAT", vdj2.sequenceFormatted)

        var vdjCombine = VdjBuilderUtils.mergeVDJs(vdj1, vdj2, VDJSequenceBuilderTest.MIN_BASE_QUALITY)

        assertNotNull(vdjCombine.vAnchor)
        assertNull(vdjCombine.jAnchor)
        assertEquals("CCTCA-GGTGATAA".replace("-", ""), vdjCombine.layout.consensusSequence())
        assertEquals("TCA-GGTGAT", vdjCombine.sequenceFormatted)
        assertEquals(3, vdjCombine.vAnchor!!.anchorBoundary)

        // try to merge it the other way, result should be the same
        vdjCombine = VdjBuilderUtils.mergeVDJs(vdj2, vdj1, VDJSequenceBuilderTest.MIN_BASE_QUALITY)

        assertNotNull(vdjCombine.vAnchor)
        assertNull(vdjCombine.jAnchor)
        assertEquals("CCTCA-GGTGATAA".replace("-", ""), vdjCombine.layout.consensusSequence())
        assertEquals("TCA-GGTGAT", vdjCombine.sequenceFormatted)
        assertEquals(3, vdjCombine.vAnchor!!.anchorBoundary)
    }

    @Test
    fun testMergeVDJsJOnly()
    {
        // layout1:        CCTCAGGTG-AT
        // VDJ1:             TCAGGTG-AT
        // layout2:            AGGTG-ATAA
        // VDJ2:               AGGTG-ATA
        // merge layout:   CCTCAGGTG-ATAA
        // merge VDJ:        TCAGGTG-ATA

        val layout1 = TestUtils.createLayout("CCTCAGGTG-AT".replace("-", ""), 5)
        val vdj1 = TestUtils.createVDJ(layout1, null, 7, 2, 11)
        assertEquals("TCAGGTG-AT", vdj1.sequenceFormatted)

        val layout2 = TestUtils.createLayout("AGGTG-ATAA".replace("-", ""), 3)
        val vdj2 = TestUtils.createVDJ(layout2, null, 5, 0, 8)
        assertEquals("AGGTG-ATA", vdj2.sequenceFormatted)

        var vdjCombine = VdjBuilderUtils.mergeVDJs(vdj1, vdj2, VDJSequenceBuilderTest.MIN_BASE_QUALITY)

        assertNull(vdjCombine.vAnchor)
        assertNotNull(vdjCombine.jAnchor)
        assertEquals("CCTCAGGTG-ATAA".replace("-", ""), vdjCombine.layout.consensusSequence())
        assertEquals("TCAGGTG-ATA", vdjCombine.sequenceFormatted)
        assertEquals(7, vdjCombine.jAnchor!!.anchorBoundary)

        // try to merge it the other way, result should be the same
        vdjCombine = VdjBuilderUtils.mergeVDJs(vdj2, vdj1, VDJSequenceBuilderTest.MIN_BASE_QUALITY)

        assertNull(vdjCombine.vAnchor)
        assertNotNull(vdjCombine.jAnchor)
        assertEquals("CCTCAGGTG-ATAA".replace("-", ""), vdjCombine.layout.consensusSequence())
        assertEquals("TCAGGTG-ATA", vdjCombine.sequenceFormatted)
        assertEquals(7, vdjCombine.jAnchor!!.anchorBoundary)
    }
}