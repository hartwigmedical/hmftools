package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.TestUtils.createReadCandidate
import com.hartwig.hmftools.cider.genes.VJ
import com.hartwig.hmftools.cider.layout.ReadLayout
import kotlin.test.*

class VdjAnnotatorTest
{
    @Test
    fun testCountNonSplitRead()
    {
        val adaptor = MockVJReadLayoutAdaptor()
        val vdjAnnotator = VdjAnnotator(adaptor)

        // length 100
        // VDJ sequence is the middle 70 bases
        // read : TGCGAATACCCACATCCTGAGAGTGTCAGA-TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACC
        // VDJ  : ------- 20 ---------GAGTGTCAGA-TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGA--- 10 ---
        //
        val readSeq = "TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACC"
        val readSliceStart = 20
        val readSliceEnd = 90
        val vdjSeq = readSeq.substring(readSliceStart, readSliceEnd)

        // V only layout
        val vdj: VDJSequence = TestUtils.createVDJ(vdjSeq, 10, null)

        val read: ReadLayout.Read = vdj.layout.reads.first()

        // set up the non split read calculation
        val readCandidate = createReadCandidate(readSeq, false,
            false, VJ.V, 0, 30)

        adaptor.readCandidateMap[read] = readCandidate
        adaptor.readSliceMap[read] = ReadSlice(readCandidate.read, false,
            readSliceStart, readSliceEnd)

        // make it non split
        //
        //
        //
        readCandidate.read.cigarString = "100M"
        assertEquals(1, vdjAnnotator.countNonSplitReads(vdj, VJ.V))

        // now try putting a split just at the start of the VDJ
        // since the read slice starts 20 bases, a split at the 20 bases
        // would not count towards non split read
        readCandidate.read.cigarString = "20S80M"
        assertEquals(1, vdjAnnotator.countNonSplitReads(vdj, VJ.V))

        // now put it inside the range
        readCandidate.read.cigarString = "21S79M"
        assertEquals(0, vdjAnnotator.countNonSplitReads(vdj, VJ.V))

        // put the split just outside the range we test for, which is 30 bases after
        // the V boundary, i.e.
        // read : TGCGAATACCCACATCCTGAGAGTGTCAGA-TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACC
        // VDJ  : ------- 20 ---------GAGTGTCAGA-TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGA--- 10 ---
        //        |-------20 ---------|-- 10---|V|-----------30---------------|
        // there for a 60M at the start would make it a non split read since it spans
        // the 10 anchor bases plus another 30 bases
        readCandidate.read.cigarString = "60M40S"
        assertEquals(1, vdjAnnotator.countNonSplitReads(vdj, VJ.V))

        // move it just 1 base less and it would not count as non split read
        readCandidate.read.cigarString = "59M41S"
        assertEquals(0, vdjAnnotator.countNonSplitReads(vdj, VJ.V))

        //
        // now test if it works for reverse complement
        //
        adaptor.readSliceMap[read] = ReadSlice(readCandidate.read, true,
            readSliceStart, readSliceEnd)

        // make it non split
        //
        //
        //
        readCandidate.read.cigarString = "100M"
        assertEquals(1, vdjAnnotator.countNonSplitReads(vdj, VJ.V))

        // now try putting a split just at the start of the VDJ
        // since the read slice starts 20 bases, a split at the 20 bases
        // would not count towards non split read
        readCandidate.read.cigarString = "80M20S"
        assertEquals(1, vdjAnnotator.countNonSplitReads(vdj, VJ.V))

        // now put it inside the range
        readCandidate.read.cigarString = "79M21S"
        assertEquals(0, vdjAnnotator.countNonSplitReads(vdj, VJ.V))

        // put the split just outside the range we test for, which is 30 bases after
        // the V boundary, i.e.
        // read : TGCGAATACCCACATCCTGAGAGTGTCAGA-TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACC
        // VDJ  : ------- 20 ---------GAGTGTCAGA-TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACCCACATCCTGAGAGTGTCAGA--- 10 ---
        //        |-------20 ---------|-- 10---|V|-----------30---------------|
        // there for a 60M at the start would make it a non split read since it spans
        // the 10 anchor bases plus another 30 bases
        readCandidate.read.cigarString = "40S60M"
        assertEquals(1, vdjAnnotator.countNonSplitReads(vdj, VJ.V))

        // move it just 1 base less and it would not count as non split read
        readCandidate.read.cigarString = "41S59M"
        assertEquals(0, vdjAnnotator.countNonSplitReads(vdj, VJ.V))
    }

    @Test
    fun testCountNonSplitReadAlignedPos()
    {
        val adaptor = MockVJReadLayoutAdaptor()
        val vdjAnnotator = VdjAnnotator(adaptor)

        // next we want to test that the reads with different aligned position is also calculated correctly
        // add another read to the layout which starts 10 bases further down
        // read1: GAGTGTCAGA-TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACC-CACATCCTGAGAGTGTCAGA
        // read2: ---------------------CACATCCTGAGAGTGTCAGATGCGAATACC-CACATCCTGAGAGTGTCAGA
        // VDJ  : GAGTGTCAGA-TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACC-CACATCCTGAGAGTGTCAGA
        val vdjSeq = "GAGTGTCAGA-TGCGAATACCCACATCCTGAGAGTGTCAGATGCGAATACC-CACATCCTGAGAGTGTCAGA".replace("-", "")

        val vdj: VDJSequence = TestUtils.createVDJ(vdjSeq, 10, 50)
        assertEquals(vdjSeq, vdj.sequence)
        val read1: ReadLayout.Read = vdj.layout.reads.first()

        val read2Seq = vdjSeq.substring(20)
        val read2: ReadLayout.Read = TestUtils.createLayoutRead("read2", read2Seq,
            "F".repeat(read2Seq.length), vdj.layout.alignedPosition - 20)
        vdj.layout.addRead(read2, 31)
        assertEquals(vdjSeq, vdj.sequence)

        // set up the non split read calculation
        val readCandidate1 = createReadCandidate(read1.sequenceString, false,
            false, VJ.V, 0, 20)
        adaptor.readCandidateMap[read1] = readCandidate1
        adaptor.readSliceMap[read1] = ReadSlice(readCandidate1.read, false,
            0, read1.sequence.size)

        val readCandidate2 = createReadCandidate(read2.sequenceString, false,
            false, VJ.V, 0, 20)

        adaptor.readCandidateMap[read2] = readCandidate2
        adaptor.readSliceMap[read2] = ReadSlice(readCandidate2.read, false,
        0, read2.sequence.size)

        // make it non split
        readCandidate2.read.cigarString = "50M"
        assertEquals(1, vdjAnnotator.countNonSplitReads(vdj, VJ.J))

        // make it a split read, since this read spans the 30 bases before J boundary to end of J anchor
        // any split within this read would make count
        readCandidate2.read.cigarString = "49M1S"
        assertEquals(0, vdjAnnotator.countNonSplitReads(vdj, VJ.J))

        readCandidate2.read.cigarString = "1S49M"
        assertEquals(0, vdjAnnotator.countNonSplitReads(vdj, VJ.J))
    }
}