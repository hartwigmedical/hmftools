package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.cdr3.layout.ReadLayout
import com.hartwig.hmftools.cdr3.layout.ReadLayoutBuilderTest
import htsjdk.samtools.SAMUtils
import org.apache.logging.log4j.Level
import org.junit.Before
import kotlin.test.*

class VDJSequenceBuilderTest
{
    @Before
    fun setUp()
    {
        org.apache.logging.log4j.core.config.Configurator.setRootLevel(Level.TRACE);
    }

    /*
    @Test
    fun testCombineVDJsSimple()
    {
        // create two VDJs, both simple read match, exactly the same
        val seq = "ATGCTGGTGT"

        val layout1 = ReadLayout()
        var seq1 = "CAGGTG"
        val baseQual1 = SAMUtils.fastqToPhred("FF::FF") // F is 37, : is 25

        // we are aligned at the T
        var read1 = ReadLayout.Read("read1", ReadKey("read1", true), seq1, baseQual1, 4)
        layout1.addRead(read1, ReadLayoutBuilderTest.MIN_BASE_QUALITY)

        val layout2 = ReadLayout()
        var seq2 = "AGGTGAT"
        val baseQual2 = SAMUtils.fastqToPhred("F:FFFF:") // F is 37, : is 25


        val vdj1 = createVDJ(seq, "1111111111", 3, 7, "vdj1-")
        val vdj2 = createVDJ(seq, "1122222211", 3, 7, "vdj2-")

        val vdjCombine = VDJSequenceBuilder.mergeVDJs(vdj1, vdj2, MIN_BASE_QUALITY)

        assertEquals(seq, vdjCombine.sequence)
        assertEquals("2233333322", vdjCombine.supportString)
    }*/

    @Test
    fun testCombineVDJsOverlap()
    {
        // want to make a VDJ seq that is  CA-GGT-GAT
        // by combining two VDJ seq, first CA-GGT-G
        // second one is                    A-GGT-GAT

        val layout1 = ReadLayout()
        val seq1 = "CAGGTG"
        val baseQual1 = SAMUtils.fastqToPhred("FF::FF") // F is 37, : is 25
        // we are aligned at the T
        val read1 = ReadLayout.Read("read1", ReadKey("read1", true), seq1, baseQual1, 4)
        layout1.addRead(read1, ReadLayoutBuilderTest.MIN_BASE_QUALITY)
        val vdj1 = createVDJ(layout1, 2, 5, 0, 6)

        val layout2 = ReadLayout()
        val seq2 = "AGGTGAT"
        val baseQual2 = SAMUtils.fastqToPhred("F:FFFF:") // F is 37, : is 25
        // aligned at the first A
        val read2 = ReadLayout.Read("read2", ReadKey("read2", true), seq2, baseQual2, 0)
        layout2.addRead(read2, ReadLayoutBuilderTest.MIN_BASE_QUALITY)
        val vdj2 = createVDJ(layout2, 1, 4, 0, 7)

        var vdjCombine = VDJSequenceBuilder.mergeVDJs(vdj1, vdj2, MIN_BASE_QUALITY)

        assertEquals("CAGGTGAT", vdjCombine.sequence)
        assertEquals(2, vdjCombine.vAnchor.anchorBoundary)
        assertEquals(5, vdjCombine.jAnchor.anchorBoundary)
        //assertEquals("1133333311", vdjCombine.supportString)

        // try to merge it the other way, result should be the same
        vdjCombine = VDJSequenceBuilder.mergeVDJs(vdj2, vdj1, MIN_BASE_QUALITY)

        assertEquals("CAGGTGAT", vdjCombine.sequence)
        assertEquals(2, vdjCombine.vAnchor.anchorBoundary)
        assertEquals(5, vdjCombine.jAnchor.anchorBoundary)
    }

    companion object
    {
        const val RADIX = 36
        const val MIN_BASE_QUALITY = 30.toByte()

        // a simple function to create VDJ sequence
        fun createVDJ(layout: ReadLayout, vAnchorBoundary: Int, jAnchorBoundary: Int,
            layoutStart: Int, layoutEnd: Int) : VDJSequence
        {
            assertTrue(layoutEnd <= layout.length)
            assertTrue(vAnchorBoundary < jAnchorBoundary)
            assertTrue(layoutStart + vAnchorBoundary < layout.length)
            assertTrue(layoutStart + jAnchorBoundary <= layout.length)

            // create VJ anchor
            val vAnchor = VJAnchorByReadMatch(
                type = VJAnchor.Type.V,
                geneType = VJGeneType.IGHV,
                anchorBoundary = vAnchorBoundary,
                matchMethod = "exact"
            )

            val jAnchor = VJAnchorByReadMatch(
                type = VJAnchor.Type.J,
                geneType = VJGeneType.IGHJ,
                anchorBoundary = jAnchorBoundary,
                matchMethod = "exact"
            )

            val vdj = VDJSequence("id", layout, layoutStart, layoutEnd, vAnchor, jAnchor)
            return vdj
        }

        /*
        // a simple function to create VDJ sequence
        fun createVDJ1(sequence: String, support: String, vAnchorBoundary: Int, jAnchorBoundary: Int, readIdPrefix: String) : VDJSequence
        {
            assertEquals(sequence.length, support.length)
            assertTrue(vAnchorBoundary < jAnchorBoundary)
            assertTrue(vAnchorBoundary < sequence.length)
            assertTrue(jAnchorBoundary <= sequence.length)

            val vdjSupport = ArrayList<VDJSequence.BaseSupport>()

            // we want to somehow intelligently create some reads to get the given support
            for (i in sequence.indices)
            {
                vdjSupport.add(VDJSequence.BaseSupport(sequence[i], createReadSet(support[i], readIdPrefix)))
            }

            // create VJ anchor
            val vAnchor = VJAnchorByReadMatch(
                type = VJAnchor.Type.V,
                geneType = VJGeneType.IGHV,
                anchorBoundary = vAnchorBoundary,
                matchMethod = "exact"
            )

            val jAnchor = VJAnchorByReadMatch(
                type = VJAnchor.Type.J,
                geneType = VJGeneType.IGHJ,
                anchorBoundary = jAnchorBoundary,
                matchMethod = "exact"
            )

            val vdj = VDJSequence("id", vdjSupport.toTypedArray(), vAnchor, jAnchor)
            return vdj
        }

        fun createReadSet(sizeCode: Char, readIdPrefix: String) : Set<ReadKey>
        {
            val readSet = HashSet<ReadKey>()

            // convert the size Code to int
            val size = sizeCode.digitToInt(radix = RADIX)
            repeat(size, { i -> readSet.add(ReadKey("${readIdPrefix}${i}", true)) })
            return readSet
        }*/
    }
}