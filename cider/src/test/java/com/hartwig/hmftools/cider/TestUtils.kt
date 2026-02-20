package com.hartwig.hmftools.cider

import com.google.common.collect.ImmutableList
import com.hartwig.hmftools.cider.genes.VJ
import com.hartwig.hmftools.cider.genes.VJAnchorTemplate
import com.hartwig.hmftools.cider.genes.VJGeneType
import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.cider.layout.TestLayoutRead
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMUtils
import java.util.concurrent.atomic.AtomicInteger
import kotlin.test.assertTrue

object TestUtils
{
    const val MIN_BASE_QUALITY = 30.toByte()

    // used to assign unique read id
    val nextReadId = AtomicInteger(1)

    fun createLayoutRead(readId: String, seq: String, baseQualityString: String, alignedPosition: Int, firstOfPair: Boolean = true)
    : ReadLayout.Read
    {
        val baseQual = SAMUtils.fastqToPhred(baseQualityString)

        // we are aligned at the T
        return TestLayoutRead(readId, ReadKey(readId, firstOfPair), seq.toByteArray(), baseQual, alignedPosition)
    }

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

        return VJReadCandidate(record, ImmutableList.of(), vjGeneType,
            "CACGTG", VJReadCandidate.MatchMethod.ALIGN,
            useReverseComplement, anchorOffsetStart, anchorOffsetEnd,
            0, 0)
    }

    // create a very simple layout with just a sequence
    fun createLayout(seq: String, alignedPosition: Int = 0, minBaseQuality: Byte = MIN_BASE_QUALITY) : ReadLayout
    {
        val baseQual = SAMUtils.phredToFastq(minBaseQuality * 2).toString().repeat(seq.length)
        val read = createLayoutRead("createLayout::autoReadId::${nextReadId.getAndIncrement()}", seq, baseQual, alignedPosition)
        val layout = ReadLayout()
        layout.addRead(read, minBaseQuality)
        return layout
    }

    fun createVDJ(layout: ReadLayout, vAnchorBoundary: Int?, jAnchorBoundary: Int?,
                  layoutStart: Int, layoutEnd: Int) : VDJSequence
    {
        assertTrue(layoutEnd <= layout.length)
        assertTrue(vAnchorBoundary == null || jAnchorBoundary == null || vAnchorBoundary < jAnchorBoundary)

        val vAnchor: VJAnchor?
        val jAnchor: VJAnchor?

        if (vAnchorBoundary != null)
        {
            assertTrue(layoutStart + vAnchorBoundary < layout.length)
            val vAnchorSeq = layout.consensusSequenceString().drop(layoutStart).substring(0, vAnchorBoundary)
            // create VJ anchor
            vAnchor = VJAnchorByReadMatch(
                vj = VJ.V,
                geneType = VJGeneType.IGHV,
                anchorBoundary = vAnchorBoundary,
                matchMethod = "exact",
                templateAnchorSeq = vAnchorSeq,
                numReads = layout.reads.size
            )
        }
        else
        {
            vAnchor = null
        }
        if (jAnchorBoundary != null)
        {
            assertTrue(layoutStart + jAnchorBoundary <= layout.length)
            val jAnchorSeq = layout.consensusSequenceString().substring(jAnchorBoundary)
            jAnchor = VJAnchorByReadMatch(
                vj = VJ.J,
                geneType = VJGeneType.IGHJ,
                anchorBoundary = jAnchorBoundary,
                matchMethod = "exact",
                templateAnchorSeq = jAnchorSeq,
                numReads = layout.reads.size
            )
        }
        else
        {
            jAnchor = null
        }

        val vdj = VDJSequence(layout, layoutStart, layoutEnd, vAnchor, jAnchor)
        return vdj
    }

    fun createVDJ(sequence: String, vAnchorBoundary: Int?, jAnchorBoundary: Int?) : VDJSequence
    {
        assertTrue(vAnchorBoundary == null || jAnchorBoundary == null || vAnchorBoundary < jAnchorBoundary)
        assertTrue(jAnchorBoundary == null || jAnchorBoundary < sequence.length)

        val layout = createLayout(sequence)
        val vAnchor: VJAnchor?
        val jAnchor: VJAnchor?

        if (vAnchorBoundary != null)
        {
            val vAnchorSeq = layout.consensusSequenceString().substring(0, vAnchorBoundary)
            // create VJ anchor
            vAnchor = VJAnchorByReadMatch(
                vj = VJ.V,
                geneType = VJGeneType.IGHV,
                anchorBoundary = vAnchorBoundary,
                matchMethod = "exact",
                templateAnchorSeq = vAnchorSeq,
                numReads = layout.reads.size
            )
        }
        else
        {
            vAnchor = null
        }
        if (jAnchorBoundary != null)
        {
            val jAnchorSeq = layout.consensusSequenceString().substring(jAnchorBoundary)
            jAnchor = VJAnchorByReadMatch(
                vj = VJ.J,
                geneType = VJGeneType.IGHJ,
                anchorBoundary = jAnchorBoundary,
                matchMethod = "exact",
                templateAnchorSeq = jAnchorSeq,
                numReads = layout.reads.size
            )
        }
        else
        {
            jAnchor = null
        }

        val vdj = VDJSequence(layout, 0, layout.length, vAnchor, jAnchor)
        return vdj
    }

    // some genes we can use for testing
    val ighJ1 = VJAnchorTemplate(
        VJGeneType.IGHJ,
        "IGHJ1",
        "01",
        null,
        "TGGGGCCAGGGCACCCTGGTCACCGTCTCC",
        null
    )

    val ighJ6 = VJAnchorTemplate(
        VJGeneType.IGHJ, "IGHJ6","01", null,
        "TGGGGGCAAGGGACCACGGTCACCGTCTCC",
        null)

    val ighV1_18 = VJAnchorTemplate(
        VJGeneType.IGHV,
        "IGHV1-18",
        "01",
        null,
        "AGATCTGACGACACGGCCGTGTATTACTGT",
        null)

    val ighV3_7 = VJAnchorTemplate(
        VJGeneType.IGHV, "IGHV3-7", "01", null,
        "AGAGCCGAGGACACGGCTGTGTATTACTGT",
        null)
}