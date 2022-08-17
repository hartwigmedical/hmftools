package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.cider.layout.ReadLayoutBuilderTest
import htsjdk.samtools.SAMUtils
import org.apache.logging.log4j.Level
import org.junit.Before
import java.util.IdentityHashMap
import kotlin.test.*

// create a mock layout adaptor that just return anchor range that we give it
class MockVJReadLayoutAdaptor : IVJReadLayoutAdaptor
{
    val anchorRangeMap = IdentityHashMap<ReadLayout, IntRange>()

    override fun getAnchorMatchType(layout: ReadLayout): VJReadCandidate.AnchorMatchMethod
    {
        return VJReadCandidate.AnchorMatchMethod.ALIGN
    }

    override fun getAnchorRange(vj: VJ, layout: ReadLayout) : IntRange?
    {
        return anchorRangeMap[layout]
    }
}

class VDJSequenceBuilderTest
{
    @Before
    fun setUp()
    {
        org.apache.logging.log4j.core.config.Configurator.setRootLevel(Level.TRACE)
    }

    // build a VDJ from a layout with just V anchor
    @Test
    fun testBuildVdjFromVLayout()
    {
        val mockVjReadLayoutAdaptor = MockVJReadLayoutAdaptor()
        val mockAnchorBlosumSearcher = MockAnchorBlosumSearcher()
        val vdjSeqBuilder = VDJSequenceBuilder(
            mockVjReadLayoutAdaptor, mockAnchorBlosumSearcher, MIN_BASE_QUALITY,
            CiderConstants.MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES)

        // we create a layout that has the V anchor but not J anchor
        // 30 bases long, V anchor is position 3-10:
        // TGC-GAATACC-CACATCCTGAGAGTGTCAGA
        //    V anchor
        // We want to set the mock anchor blosum searcher to return J anchor at position 20-25:
        // TGC-GAATACC-CACATCCTGA-GAGTG-TCAGA
        //     |_____|            |___|
        //    V anchor(align)    J anchor(blosum)

        val seq = "TGC-GAATACC-CACATCCTGA-GAGTG-TCAGA".replace("-", "")
        val layout = TestUtils.createLayout(seq)

        mockVjReadLayoutAdaptor.anchorRangeMap[layout] = 3 until 10

        // create a blosum match of type IGHJ, which will match with IGHV
        mockAnchorBlosumSearcher.anchorBlosumMatch = AnchorBlosumMatch(20, 25, "CCCCC",
            listOf(TestUtils.ighJ1), 10)

        // now we should be able to build a VDJ sequence
        val vdjSeq = vdjSeqBuilder.tryCompleteLayoutWithBlosum(VJGeneType.IGHV, layout, 0)

        assertNotNull(vdjSeq)
        assertEquals(VJGeneType.IGHV, vdjSeq.vAnchor.geneType)
        assertEquals(VJGeneType.IGHJ, vdjSeq.jAnchor.geneType)
        assertIs<VJAnchorByReadMatch>(vdjSeq.vAnchor)
        assertIs<VJAnchorByBlosum>(vdjSeq.jAnchor)
        assertEquals("GAATACC", vdjSeq.vAnchorSequence)
        assertEquals("GAGTG", vdjSeq.jAnchorSequence)
        assertEquals("CACATCCTGA", vdjSeq.cdr3SequenceShort)
    }

    // build a VDJ from a layout with just J anchor
    @Test
    fun testBuildVdjFromJLayout()
    {
        val mockVjReadLayoutAdaptor = MockVJReadLayoutAdaptor()
        val mockAnchorBlosumSearcher = MockAnchorBlosumSearcher()
        val vdjSeqBuilder = VDJSequenceBuilder(
            mockVjReadLayoutAdaptor, mockAnchorBlosumSearcher, MIN_BASE_QUALITY,
            CiderConstants.MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES)

        // we create a layout that has the J anchor but not V anchor
        // 30 bases long, J anchor is position 20-25:
        // We want to set the mock anchor blosum searcher to return V anchor at position 3-10:
        // TGC-GAATACC-CACATCCTGA-GAGTG-TCAGA
        //     |_____|            |___|
        //    V anchor(blosum)    J anchor(align)

        val seq = "TGC-GAATACC-CACATCCTGA-GAGTG-TCAGA".replace("-", "")
        val layout = TestUtils.createLayout(seq)

        mockVjReadLayoutAdaptor.anchorRangeMap[layout] = 20 until 25

        // create a blosum match of type IGHV, which will match with IGHJ
        mockAnchorBlosumSearcher.anchorBlosumMatch = AnchorBlosumMatch(3, 10, "CCCCC",
            listOf(TestUtils.ighV1_18), 10)

        // now we should be able to build a VDJ sequence
        val vdjSeq = vdjSeqBuilder.tryCompleteLayoutWithBlosum(VJGeneType.IGHJ, layout, 0)

        assertNotNull(vdjSeq)
        assertEquals(VJGeneType.IGHV, vdjSeq.vAnchor.geneType)
        assertEquals(VJGeneType.IGHJ, vdjSeq.jAnchor.geneType)
        assertIs<VJAnchorByBlosum>(vdjSeq.vAnchor)
        assertIs<VJAnchorByReadMatch>(vdjSeq.jAnchor)
        assertEquals("GAATACC", vdjSeq.vAnchorSequence)
        assertEquals("GAGTG", vdjSeq.jAnchorSequence)
        assertEquals("CACATCCTGA", vdjSeq.cdr3SequenceShort)
    }

    // test building VDJ sequence by merging overlapping V layout and J layout
    @Test
    fun testBuildVdjFromOverlapLayouts()
    {
        val mockVjReadLayoutAdaptor = MockVJReadLayoutAdaptor()
        val mockAnchorBlosumSearcher = MockAnchorBlosumSearcher()
        val vdjSeqBuilder = VDJSequenceBuilder(
            mockVjReadLayoutAdaptor, mockAnchorBlosumSearcher, MIN_BASE_QUALITY,
            10)

        // we create two layouts, a V layout and a J layout, and they overlap each other by 13 bases:
        // v layout:   TGC-GAATACC-CACATCCTGA-G
        // j layout:            CC-CACATCCTGA-GAGTG-TCAGA
        //                 |_____|            |___|
        //            V anchor(align)    J anchor(align)
        //
        // V anchor is position 3-10 in the V layout, J anchor is position 12-17 in the J layout.
        // The final VDJ is 30 bases long, V anchor at position 3-10, J anchor is position 20-25

        val seq = "TGC-GAATACC-CACATCCTGA-GAGTG-TCAGA".replace("-", "")
        val vLayout = TestUtils.createLayout(seq.substring(0, 21))
        val jLayout = TestUtils.createLayout(seq.substring(8))

        mockVjReadLayoutAdaptor.anchorRangeMap[vLayout] = 3 until 10
        mockVjReadLayoutAdaptor.anchorRangeMap[jLayout] = 12 until 17

        // now we should be able to build a VDJ sequence
        val vdjSeq = vdjSeqBuilder.tryOverlapVJ(vLayout, jLayout, VJGeneType.TRAV, VJGeneType.TRAJ)

        assertNotNull(vdjSeq)
        assertEquals(VJGeneType.TRAV, vdjSeq.vAnchor.geneType)
        assertEquals(VJGeneType.TRAJ, vdjSeq.jAnchor.geneType)
        assertIs<VJAnchorByReadMatch>(vdjSeq.vAnchor)
        assertIs<VJAnchorByReadMatch>(vdjSeq.jAnchor)
        assertEquals("GAATACC", vdjSeq.vAnchorSequence)
        assertEquals("GAGTG", vdjSeq.jAnchorSequence)
        assertEquals("CACATCCTGA", vdjSeq.cdr3SequenceShort)
    }

    @Test
    fun testMergeVDJsSimple()
    {
        // create two VDJs, both simple read match, exactly the same
        val seq = "ATGCTGGTGT"

        val layout1 = ReadLayout()

        // we are aligned at the T
        val baseQual1 = SAMUtils.fastqToPhred("FF:FFFF:FF") // F is 37, : is 25
        val read1 = ReadLayout.Read("read1", ReadKey("read1", true), seq, baseQual1, 4)
        layout1.addRead(read1, ReadLayoutBuilderTest.MIN_BASE_QUALITY)

        val layout2 = ReadLayout()
        val baseQual2 = SAMUtils.fastqToPhred("FFFF::FFFF") // F is 37, : is 25
        val read2 = ReadLayout.Read("read2", ReadKey("read2", true), seq, baseQual2, 4)
        layout2.addRead(read2, ReadLayoutBuilderTest.MIN_BASE_QUALITY)

        val vdj1 = createVDJ(layout1, 3, 7, 0, 10)
        val vdj2 = createVDJ(layout2, 3, 7, 0, 10)

        val vdjCombine = VDJSequenceBuilder.mergeVDJs(vdj1, vdj2, MIN_BASE_QUALITY)

        assertEquals(seq, vdjCombine.sequence)
        assertEquals("2212112122", vdjCombine.supportString) // due to low base qual in some reads
    }

    @Test
    fun testMergeVDJsLongShort()
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