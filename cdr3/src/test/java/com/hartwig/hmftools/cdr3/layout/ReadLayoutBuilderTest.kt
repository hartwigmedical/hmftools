package com.hartwig.hmftools.cdr3.layout

import com.hartwig.hmftools.cdr3.ReadKey
import htsjdk.samtools.SAMUtils
import kotlin.test.*

class ReadLayoutBuilderTest
{
    companion object
    {
        const val MIN_BASE_QUALITY = 30.toByte()
        const val MIN_OVERLAP_BASES = 3
    }

    @Test
    fun testCompareReads1()
    {
        val seq1 = "CAGGTGCAGCTGGTGGAGTCTGGGGGA"
        val seq2 = "GAGGTGCAGCTGGTAGAGTCTGGGAGA"
        val baseQual1 = ByteArray(seq1.length){ 50 }
        val baseQual2 = baseQual1

        val (matchCount, compareCount) = ReadLayoutBuilder.sequenceMatchCount(seq1, seq2, baseQual1, baseQual2,
            0, 0, 20, 30)

        assertEquals(18, matchCount)
    }

    @Test
    fun testCompareReads2()
    {
        val seq1 = "CAGGTGCAGCTGGTGGAGTCTGGGGGA"
        val seq2 = "GTGCAGCTGGTAGAGTCTGGGAGA"
        val baseQual1 = ByteArray(seq1.length){ 50 }
        val baseQual2 = baseQual1

        val (matchCount, compareCount) = ReadLayoutBuilder.sequenceMatchCount(seq1, seq2, baseQual1, baseQual2,
            3, 0, 20, 30)

        assertEquals(19, matchCount)
    }

    @Test
    fun testAddToOverlay()
    {
        val group = ReadLayout()
        var seq = "CAGGTG"

        // we are aligned at the T
        var readData = ReadLayout.Read("read1", ReadKey("read1", true), seq, ByteArray(seq.length){ 50 }, 4)
        ReadLayoutBuilder.addToOverlay(group, readData, 30)

        assertEquals(seq, group.consensusSequence())
        assertEquals(4, group.alignedPosition)

        // now test a sequence that extend the start by 3 bases
        seq = "AGCCAGGT"
        readData = ReadLayout.Read("read2", ReadKey("read2", true), seq, ByteArray(seq.length){ 50 }, 7)
        ReadLayoutBuilder.addToOverlay(group, readData, 30)

        assertEquals(seq + "G", group.consensusSequence())
        assertEquals(7, group.alignedPosition)

        // now test a sequence that extend the end by 3 bases
        seq = "AGGTGCAA"
        readData = ReadLayout.Read("read3", ReadKey("read3", true), seq, ByteArray(seq.length){ 50 }, 3)
        ReadLayoutBuilder.addToOverlay(group, readData, 30)

        assertEquals("AGCC" + seq, group.consensusSequence())
        assertEquals(7, group.alignedPosition)

        assertEquals(3, group.reads.size)

        // lets see if we got the correct support counts
        assertEquals(mapOf('A' to 1), group.highQualSequenceSupport.support[0].countMap)
        assertEquals(mapOf('G' to 1), group.highQualSequenceSupport.support[1].countMap)
        assertEquals(mapOf('C' to 1), group.highQualSequenceSupport.support[2].countMap)
        assertEquals(mapOf('C' to 2), group.highQualSequenceSupport.support[3].countMap)
        assertEquals(mapOf('A' to 3), group.highQualSequenceSupport.support[4].countMap)
        assertEquals(mapOf('G' to 3), group.highQualSequenceSupport.support[5].countMap)
        assertEquals(mapOf('G' to 3), group.highQualSequenceSupport.support[6].countMap)
        assertEquals(mapOf('T' to 3), group.highQualSequenceSupport.support[7].countMap)
        assertEquals(mapOf('G' to 2), group.highQualSequenceSupport.support[8].countMap)
        assertEquals(mapOf('C' to 1), group.highQualSequenceSupport.support[9].countMap)
        assertEquals(mapOf('A' to 1), group.highQualSequenceSupport.support[10].countMap)
        assertEquals(mapOf('A' to 1), group.highQualSequenceSupport.support[11].countMap)
    }

    @Test
    fun testLayoutMatch()
    {
        val group = ReadLayout()
        var seq = "CAGGTG"

        // we are aligned at the T
        var readData = ReadLayout.Read("read1", ReadKey("read1", true), seq, ByteArray(seq.length){ 50 }, 4)
        ReadLayoutBuilder.addToOverlay(group, readData, 30)

        assertEquals(seq, group.consensusSequence())
        assertEquals(4, group.alignedPosition)

        // match a new sequence against the overlay
        seq = "AGCCAGAT"
        readData = ReadLayout.Read("read2", ReadKey("read2", true), seq, ByteArray(seq.length){ 50 }, 7)

        //val (matchCount, compareCount) = ReadLayoutBuilder.layoutMatchCount(group, readData, false, 30)

        //assertEquals(4, matchCount)
        //assertEquals(5, compareCount)
    }

    @Test
    fun testLayoutSupport()
    {
        val group = ReadLayout()
        var seq = "CAGGTG"

        // we are aligned at the T
        var readData = ReadLayout.Read("read1", ReadKey("read1", true), seq, ByteArray(seq.length){ 50 }, 4)
        ReadLayoutBuilder.addToOverlay(group, readData, 30)

        assertEquals(seq, group.consensusSequence())
        assertEquals(4, group.alignedPosition)

        // match a new sequence against the overlay
        seq = "AGCCAGAT"
        readData = ReadLayout.Read("read2", ReadKey("read2", true), seq, ByteArray(seq.length){ 50 }, 7)

        //val (matchCount, compareCount) = ReadLayoutBuilder.layoutMatchCount(group, readData, false, 30)

        //assertEquals(4, matchCount)
        //assertEquals(5, compareCount)
    }

    @Test
    fun testLayoutMerge()
    {
        val layout1 = ReadLayout()
        var seq1 = "CAGGTG"
        val baseQual1 = SAMUtils.fastqToPhred("FF::FF") // F is 37, : is 25

        // we are aligned at the T
        var read1 = ReadLayout.Read("read1", ReadKey("read1", true), seq1, baseQual1, 4)
        layout1.addRead(read1, MIN_BASE_QUALITY)

        val layout2 = ReadLayout()
        var seq2 = "AGGTGAT"
        val baseQual2 = SAMUtils.fastqToPhred("F:FFFF:") // F is 37, : is 25

        // aligned at the first A
        var read2 = ReadLayout.Read("read2", ReadKey("read2", true), seq2, baseQual2, 0)
        layout2.addRead(read2, MIN_BASE_QUALITY)

        // now we have 2 layouts, one is GAGGTG, another is AGGTGAT, we merge them together to create CAGGTGAT
        // and we want to align them at the 2nd G
        val layout3 = ReadLayout.merge(layout1, layout2, -1, 2, MIN_BASE_QUALITY)

        assertEquals("CAGGTGAT", layout3.consensusSequence())
        assertEquals("12012210", layout3.highQualSupportString())
        //assertEquals(5, compareCount)
    }
}
