package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.cdr3.layout.ReadLayout
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter
import htsjdk.samtools.SAMUtils
import java.io.BufferedWriter

object VJReadOverlayFile
{
    private const val EXTENSION = ".overlay"

    @JvmStatic
    fun writeOverlays(filePath: String, overlayMap: Map<VJGeneType, List<ReadLayout>>, minSoftClipBases: Int)
    {
        createBufferedWriter(filePath).use { writer ->
            for ((geneType, overlayList) in overlayMap)
            {
                writeOverlays(writer, geneType, overlayList, minSoftClipBases)
            }
        }
    }

    private fun writeOverlays(writer: BufferedWriter, geneType: VJGeneType, overlayList: List<ReadLayout>, minSoftClipBases: Int)
    {
        // sort the overlays by number of reads
        val sortedOverlayList = overlayList.sortedByDescending({ o -> o.reads.size })

        for (i in sortedOverlayList.indices)
        {
            writeLayout1(sortedOverlayList, i, geneType, writer)
        }
    }

    private fun writeLayout(
        sortedOverlayList: List<ReadLayout>,
        i: Int,
        geneType: VJGeneType,
        writer: BufferedWriter
    )
    {
        val overlay: ReadLayout = sortedOverlayList[i]

        // count how many split reads
        val numSplitReads5Bases = overlay.reads.map { o -> VJReadLayoutAdaptor.toReadCandidate(o) }
            .count { o -> Math.max(o.leftSoftClip, o.rightSoftClip) >= 5 }

        val numSplitReads10Bases = overlay.reads.map { o -> VJReadLayoutAdaptor.toReadCandidate(o) }
            .count { o -> Math.max(o.leftSoftClip, o.rightSoftClip) >= 10 }

        val anchorRange = VJReadLayoutAdaptor.getAnchorRange(geneType, overlay)!!
        val cdr3Range = VJReadLayoutAdaptor.getCdr3Range(geneType, overlay)!!

        val range1 = if (geneType.isV) anchorRange else cdr3Range
        val range2 = if (geneType.isV) cdr3Range else anchorRange

        val sequence1 = safeSubstring(overlay.consensusSequence(), range1)
        val sequence2 = safeSubstring(overlay.consensusSequence(), range2)

        val support1 = safeSubstring(overlay.highQualSupportString(), range1)
        val support2 = safeSubstring(overlay.highQualSupportString(), range2)

        // make sure aligned to codon
        val aa1 = Codons.aminoAcidFromBases(sequence1.substring(sequence1.length % 3))
        val aa2 = Codons.aminoAcidFromBases(sequence2)

        writer.write("${i} type: ${geneType}, read count: ${overlay.reads.size}, split(5) read count: ${numSplitReads5Bases}, ")
        writer.write("split(10) read count: ${numSplitReads10Bases}, ")
        writer.write("AA: ${aa1}-${aa2}\n")
        writer.write("    ${sequence1}-${sequence2}\n")
        writer.write("    ${support1}-${support2}\n")

        for (r in overlay.reads)
        {
            val read: VJReadCandidate = VJReadLayoutAdaptor.toReadCandidate(r)
            val readPadding = Math.max(overlay.alignedPosition - r.alignedPosition, 0)
            val paddedSeq = " ".repeat(readPadding) + read.readSequence
            val paddedQual = " ".repeat(readPadding) + read.baseQualityString
            writer.write("    read: ${read.read}\n")
            writer.write("    ${safeSubstring(paddedSeq, range1)}-${safeSubstring(paddedSeq, range2)}\n")
            writer.write("    ${safeSubstring(paddedQual, range1)}-${safeSubstring(paddedQual, range2)}\n")
        }
    }

    private fun writeLayout1(
        sortedOverlayList: List<ReadLayout>,
        i: Int,
        geneType: VJGeneType,
        writer: BufferedWriter
    )
    {
        val overlay: ReadLayout = sortedOverlayList[i]

        // count how many split reads
        val numSplitReads5Bases = overlay.reads.map { o -> VJReadLayoutAdaptor.toReadCandidate(o) }
            .count { o -> Math.max(o.leftSoftClip, o.rightSoftClip) >= 5 }

        val numSplitReads10Bases = overlay.reads.map { o -> VJReadLayoutAdaptor.toReadCandidate(o) }
            .count { o -> Math.max(o.leftSoftClip, o.rightSoftClip) >= 10 }

        val anchorRange = VJReadLayoutAdaptor.getAnchorRange(geneType, overlay)!!
        var sequence = overlay.consensusSequence()
        var support = overlay.highQualSupportString()

        // make sure aligned to codon
        val codonAlignedSeq = sequence.drop(anchorRange.first % 3)
        val aa = insertDashes(Codons.aminoAcidFromBases(codonAlignedSeq), IntRange(anchorRange.first / 3, (anchorRange.last + 1) / 3 - 1))

        // insert - into the sequence and support
        sequence = insertDashes(sequence, anchorRange)
        support = insertDashes(support, anchorRange)

        writer.write("${i} type: ${geneType}, read count: ${overlay.reads.size}, split(5) read count: ${numSplitReads5Bases}, ")
        writer.write("split(10) read count: ${numSplitReads10Bases}, ")
        writer.write("AA: ${aa}\n")
        writer.write("    ${sequence}\n")
        writer.write("    ${support}\n")

        for (r in overlay.reads)
        {
            val read: VJReadCandidate = VJReadLayoutAdaptor.toReadCandidate(r)
            val readPadding = Math.max(overlay.alignedPosition - r.alignedPosition, 0)
            val paddedSeq = " ".repeat(readPadding) + r.sequence
            val paddedQual = " ".repeat(readPadding) + SAMUtils.phredToFastq(r.baseQualities)
            writer.write("    read: ${read.read}\n")
            writer.write("    ${insertDashes(paddedSeq, anchorRange)}\n")
            writer.write("    ${insertDashes(paddedQual, anchorRange)}\n")
        }
    }

    private fun safeSubstring(str: String, range: IntRange): String
    {
        if (str.length <= range.first)
        {
            return ""
        }
        if (str.length <= range.last)
        {
            return str.substring(range.first)
        }
        return str.substring(range)
    }

    private fun insertDashes(str: String, anchorRange: IntRange): String
    {
        if (str.length <= anchorRange.first)
        {
            return str
        }
        if (str.length <= anchorRange.last)
        {
            return str.substring(0, anchorRange.first) + '-' + str.substring(anchorRange.first)
        }
        if (anchorRange.first == 0)
        {
            // omit the dash if it is at the start
            return str.substring(anchorRange) + '-' + str.substring(anchorRange.last + 1)
        }
        if (anchorRange.last + 1 == str.length)
        {
            // omit the dash at the end
            return str.substring(0, anchorRange.first) + '-' + str.substring(anchorRange)
        }
        return str.substring(0, anchorRange.first) + '-' + str.substring(anchorRange) + '-' + str.substring(anchorRange.last + 1)
    }
}