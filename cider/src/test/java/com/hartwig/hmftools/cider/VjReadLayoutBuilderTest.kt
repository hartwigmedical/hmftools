package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.TestUtils.createReadCandidate
import com.hartwig.hmftools.cider.genes.VJ
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
        assertEquals("AAGACACGGC", layoutRead.sequenceString)
        assertEquals(5, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5
        anchorOffsetEnd = 8
        readCandidate = createReadCandidate(seq, false,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequenceString)
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
        assertEquals("AAGACACGGC", layoutRead.sequenceString)
        assertEquals(5, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5
        anchorOffsetEnd = 8
        readCandidate = createReadCandidate(seq, false,
            true, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequenceString)
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
        assertEquals("AAGACACGGC", layoutRead.sequenceString)
        assertEquals(5, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5
        anchorOffsetEnd = 8
        readCandidate = createReadCandidate(seq, false,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequenceString)
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
        assertEquals("AAGACACGGC", layoutRead.sequenceString) // remove the poly C plus another 5 bases
        assertEquals(5, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5 + 13
        anchorOffsetEnd = 8 + 13
        readCandidate = createReadCandidate(seq, true,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequenceString) // remove poly C plus another 5 bases
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
        assertEquals("AAGACACGGC", layoutRead.sequenceString) // remove poly G plus another 5 bases
        assertEquals(5, layoutRead.alignedPosition) // aligned at last base of anchor

        // now test J read
        anchorOffsetStart = 5 + 1
        anchorOffsetEnd = 8 + 1
        readCandidate = createReadCandidate(seq, false,
            false, VJ.J, anchorOffsetStart, anchorOffsetEnd)

        layoutRead = layoutAdaptor.readCandidateToLayoutRead(readCandidate)
        assertNotNull(layoutRead)
        assertEquals("AAGACACGGC", layoutRead.sequenceString) // remove ploy G plus another 5 bases
        assertEquals(5, layoutRead.alignedPosition) // aligned at first base of anchor
    }
}