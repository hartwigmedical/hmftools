package com.hartwig.hmftools.cider

import htsjdk.samtools.SAMRecord
import org.eclipse.collections.api.factory.Lists
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotNull

class VjReadLayoutAdaptorTest
{
    @Test
    fun testReadCandidateToLayoutRead()
    {
        val layoutAdaptor = VJReadLayoutAdaptor(0)
        val seq = "AAGACACGGC" // 10 bases long

        // test V read
        var anchorOffsetStart = 2
        var anchorOffsetEnd = 6
        var readCandidate = createReadCandidate(seq, false,
            false, VJ.V, anchorOffsetStart, anchorOffsetEnd)

        var layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)

        assertNotNull(layoutRead)
        assertEquals("GACACGGC", layoutRead.sequence) // remove the bases before anchor
        assertEquals(3, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5
        anchorOffsetEnd = 8
        readCandidate = createReadCandidate(seq, false,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACG", layoutRead.sequence) // remove the bases after anchor
        assertEquals(5, layoutRead.alignedPosition) // aligned at first base of anchor
    }

    @Test
    fun testReadCandidateToLayoutReadPolyG()
    {
        val layoutAdaptor = VJReadLayoutAdaptor(0)
        val seq = "AAGACACGGC" + "GTACT" + "G".repeat(8) // 10 bases long + 5 bases + 8 Gs

        // test V read
        var anchorOffsetStart = 2
        var anchorOffsetEnd = 6
        var readCandidate = createReadCandidate(seq, false,
            false, VJ.V, anchorOffsetStart, anchorOffsetEnd)

        var layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)

        assertNotNull(layoutRead)
        assertEquals("GACACGGC", layoutRead.sequence) // remove the bases before anchor
        assertEquals(3, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5
        anchorOffsetEnd = 8
        readCandidate = createReadCandidate(seq, false,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACG", layoutRead.sequence) // remove the bases after anchor
        assertEquals(5, layoutRead.alignedPosition) // aligned at first base of anchor
    }

    @Test
    fun testReadCandidateToLayoutReadPolyC()
    {
        val layoutAdaptor = VJReadLayoutAdaptor(0)
        val seq = "C".repeat(8) + "GTACT" + "AAGACACGGC" //  8 Cs + 5 bases + 10 bases long

        // test V read
        var anchorOffsetStart = 2 + 13 // 8 Cs + 5 bases
        var anchorOffsetEnd = 6 + 13
        var readCandidate = createReadCandidate(seq, true,
            false, VJ.V, anchorOffsetStart, anchorOffsetEnd)

        var layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)

        assertNotNull(layoutRead)
        assertEquals("GACACGGC", layoutRead.sequence) // remove the bases before anchor
        assertEquals(3, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5 + 13
        anchorOffsetEnd = 8 + 13
        readCandidate = createReadCandidate(seq, true,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACG", layoutRead.sequence) // remove the bases after anchor
        assertEquals(5, layoutRead.alignedPosition) // aligned at first base of anchor
    }

    @Test
    fun testReadCandidateToLayoutReadTrimBasesPolyG()
    {
        val layoutAdaptor = VJReadLayoutAdaptor(1)
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
        assertEquals("GACACGGC", layoutRead.sequence) // remove the bases before anchor
        assertEquals(3, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5 + 1
        anchorOffsetEnd = 8 + 1
        readCandidate = createReadCandidate(seq, false,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACG", layoutRead.sequence) // remove the bases after anchor
        assertEquals(5, layoutRead.alignedPosition) // aligned at first base of anchor
    }

    @Test
    fun testReadCandidateToLayoutReadTrimBasesPolyC()
    {

    }

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
            "CACGTG", VJReadCandidate.AnchorMatchMethod.ALIGN,
                                useReverseComplement, anchorOffsetStart, anchorOffsetEnd,
                                null, 0, 0)
    }
}