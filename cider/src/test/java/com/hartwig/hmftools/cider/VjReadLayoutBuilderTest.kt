package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.TestUtils.MIN_BASE_QUALITY
import com.hartwig.hmftools.cider.layout.TestLayoutRead
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMUtils
import org.eclipse.collections.api.factory.Lists
import kotlin.test.*

class VjReadLayoutBuilderTest
{
    @Test
    fun testReadCandidateToLayoutReadPosStrand()
    {
        val layoutAdaptor = VJReadLayoutBuilder(0, 30)
        val seq = "AAGACACGGC" // 10 bases long

        // test V read
        var anchorOffsetStart = 2
        var anchorOffsetEnd = 6
        var readCandidate = createReadCandidate(seq, false,
            false, VJ.V, anchorOffsetStart, anchorOffsetEnd)

        var layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)

        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequence)
        assertEquals(5, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5
        anchorOffsetEnd = 8
        readCandidate = createReadCandidate(seq, false,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequence)
        assertEquals(5, layoutRead.alignedPosition) // aligned at first base of anchor
    }

    @Test
    fun testReadCandidateToLayoutReadNegStrand()
    {
        val layoutAdaptor = VJReadLayoutBuilder(0, 30)
        val seq = "GCCGTGTCTT" // 10 bases long, reverse comp of AAGACACGGC

        // test V read
        var anchorOffsetStart = 2
        var anchorOffsetEnd = 6
        var readCandidate = createReadCandidate(seq, false,
            true, VJ.V, anchorOffsetStart, anchorOffsetEnd)

        var layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)

        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequence)
        assertEquals(5, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5
        anchorOffsetEnd = 8
        readCandidate = createReadCandidate(seq, false,
            true, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequence)
        assertEquals(5, layoutRead.alignedPosition) // aligned at first base of anchor
    }

    @Test
    fun testReadCandidateToLayoutReadPolyG()
    {
        val layoutAdaptor = VJReadLayoutBuilder(0, 30)
        val seq = "AAGACACGGC" + "GTACT" + "G".repeat(8) // 10 bases long + 5 bases + 8 Gs

        // test V read
        var anchorOffsetStart = 2
        var anchorOffsetEnd = 6
        var readCandidate = createReadCandidate(seq, false,
            false, VJ.V, anchorOffsetStart, anchorOffsetEnd)

        var layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)

        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequence)
        assertEquals(5, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5
        anchorOffsetEnd = 8
        readCandidate = createReadCandidate(seq, false,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequence)
        assertEquals(5, layoutRead.alignedPosition) // aligned at first base of anchor
    }

    @Test
    fun testReadCandidateToLayoutReadPolyC()
    {
        val layoutAdaptor = VJReadLayoutBuilder(0, 30)
        val seq = "C".repeat(8) + "GTACT" + "AAGACACGGC" //  8 Cs + 5 bases + 10 bases long

        // test V read
        var anchorOffsetStart = 2 + 13 // 8 Cs + 5 bases
        var anchorOffsetEnd = 6 + 13
        var readCandidate = createReadCandidate(seq, true,
            false, VJ.V, anchorOffsetStart, anchorOffsetEnd)

        var layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)

        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequence) // remove the poly C plus another 5 bases
        assertEquals(5, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5 + 13
        anchorOffsetEnd = 8 + 13
        readCandidate = createReadCandidate(seq, true,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequence) // remove poly C plus another 5 bases
        assertEquals(5, layoutRead.alignedPosition) // aligned at first base of anchor
    }

    @Test
    fun testReadCandidateToLayoutReadTrimBasesPolyG()
    {
        val layoutAdaptor = VJReadLayoutBuilder(1, 30)
        var seq = "AAGACACGGC" + "GTACT" + "G".repeat(8) // 10 bases long + 5 bases + 8 Gs

        // add 1 base around for trim
        seq = "T" + seq + "A"

        // test V read
        var anchorOffsetStart = 2 + 1
        var anchorOffsetEnd = 6 + 1
        var readCandidate = createReadCandidate(seq, false,
            false, VJ.V, anchorOffsetStart, anchorOffsetEnd)

        var layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)

        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequence) // remove poly G plus another 5 bases
        assertEquals(5, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5 + 1
        anchorOffsetEnd = 8 + 1
        readCandidate = createReadCandidate(seq, false,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequence) // remove ploy G plus another 5 bases
        assertEquals(5, layoutRead.alignedPosition) // aligned at first base of anchor
    }

    // test the function that determine if the read belongs to this layout
    @Test
    fun testReadMatchesLayout()
    {
        val layoutSeq = "TGCGAATACCCACATCCTGAGAGTGTCAGA"
        val layout = TestUtils.createLayout(layoutSeq)
        layout.alignedPosition = 15

        val readSeq1 = StringBuilder(layoutSeq.drop(4))
        val baseQual1 = SAMUtils.fastqToPhred("F".repeat(readSeq1.length)) // F is 37, : is 25
        var read = TestLayoutRead("read1", ReadKey("read1", true), readSeq1.toString(), baseQual1, 11)

        // they should match up
        assertTrue(VJReadLayoutBuilder.readMatchesLayout(read, layout, 5, MIN_BASE_QUALITY,
            10, 10))

        // test that negative value correctly disable the feature
        assertTrue(VJReadLayoutBuilder.readMatchesLayout(read, layout, 5, MIN_BASE_QUALITY,
            -1, -1))

        // change one base before aligned pos
        val readSeq2 = StringBuilder(readSeq1.toString())
        readSeq2[4] = 'T'
        read = TestLayoutRead("read1", ReadKey("read1", true), readSeq2.toString(), baseQual1, 11)

        // now cannot match
        assertFalse(VJReadLayoutBuilder.readMatchesLayout(read, layout, 5, MIN_BASE_QUALITY,
            10, 10))

        // also cannot match if we use -1 to compare everything before aligned position
        assertFalse(VJReadLayoutBuilder.readMatchesLayout(read, layout, 5, MIN_BASE_QUALITY,
            -1, 10))

        // but we can if we extend the number of bases to not compare
        // aligned pos is 11, the index 4 is changed, so we just need to only compare 6 bases
        //
        // 0123456789AB
        //      ------  <--- can only compare 6 bases
        assertTrue(VJReadLayoutBuilder.readMatchesLayout(read, layout, 5, MIN_BASE_QUALITY,
            6, 10))

        // change one base after aligned pos
        val readSeq3 = StringBuilder(readSeq1.toString())
        readSeq3[17] = 'T'
        read = TestLayoutRead("read1", ReadKey("read1", true), readSeq3.toString(), baseQual1, 11)

        // now cannot match
        assertFalse(VJReadLayoutBuilder.readMatchesLayout(read, layout, 5, MIN_BASE_QUALITY,
            10, 10))

        // also cannot match if we use -1 to compare everything after aligned position
        assertFalse(VJReadLayoutBuilder.readMatchesLayout(read, layout, 5, MIN_BASE_QUALITY,
            10, -1))

        // but we can if we extend the number of bases to not compare
        // aligned pos is 11, the index 17 is changed, so we just need to only compare 6 bases
        //
        // 111111111111
        // 123456789
        // ------  <--- can only compare 6 bases
        assertTrue(VJReadLayoutBuilder.readMatchesLayout(read, layout, 5, MIN_BASE_QUALITY,
            10, 6))

        // test that low qual mismatch is allowed
        val readSeq4 = StringBuilder(readSeq1.toString())
        readSeq4[4] = 'T'
        readSeq4[18] = 'T'
        baseQual1[4] = (MIN_BASE_QUALITY - 1).toByte()
        baseQual1[18] = (MIN_BASE_QUALITY - 1).toByte()
        read = TestLayoutRead("read1", ReadKey("read1", true), readSeq4.toString(), baseQual1, 11)

        // still match
        assertTrue(VJReadLayoutBuilder.readMatchesLayout(read, layout, 5, MIN_BASE_QUALITY,
            -1, -1))

        // but will not match if I lower the min base quality
        assertFalse(VJReadLayoutBuilder.readMatchesLayout(read, layout, 5, (MIN_BASE_QUALITY - 1).toByte(),
            -1, -1))
    }

    companion object
    {
        // helper function to create read candidates
        fun createReadCandidate(seq: String, isReadNegativeStrand: Boolean, useReverseComplement: Boolean, vj: VJ,
                                anchorOffsetStart: Int, anchorOffsetEnd: Int) : VJReadCandidate
        {
            val record = SAMRecord(null)
            record.readName = "read"
            record.readPairedFlag = true
            record.firstOfPairFlag = true
            record.readNegativeStrandFlag = isReadNegativeStrand
            record.readString = seq
            record.baseQualityString = "F".repeat(seq.length)
            val vjGeneType = if (vj == VJ.V) VJGeneType.TRAV else VJGeneType.TRAJ

            return VJReadCandidate(record, Lists.immutable.empty(), vjGeneType,
                "CACGTG", VJReadCandidate.MatchMethod.ALIGN,
                                    useReverseComplement, anchorOffsetStart, anchorOffsetEnd,
                                    0, 0)
        }
    }
}