package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.cdr3.layout.ReadLayout
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter
import htsjdk.samtools.SAMUtils
import java.io.BufferedWriter
import java.io.File

object VJReadLayoutFile
{
    private const val FILE_EXTENSION = ".cider.layout"

    private fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun writeLayouts(outputDir: String, sampleId: String, overlayMap: Map<VJGeneType, List<ReadLayout>>, minSoftClipBases: Int)
    {
        val filePath = generateFilename(outputDir, sampleId)

        createBufferedWriter(filePath).use { writer ->
            for ((geneType, overlayList) in overlayMap)
            {
                writeLayouts(writer, geneType, overlayList)
            }
        }
    }

    private fun writeLayouts(writer: BufferedWriter, geneType: VJGeneType, overlayList: List<ReadLayout>)
    {
        val vjReadLayoutAdaptor = VJReadLayoutAdaptor()

        // sort the overlays by number of reads
        val sortedOverlayList = overlayList.sortedByDescending({ o -> o.reads.size })

        for (layout in sortedOverlayList)
        {
            writeLayout(layout, geneType, vjReadLayoutAdaptor, writer)
        }
    }

    private fun writeLayout(
        overlay: ReadLayout,
        geneType: VJGeneType,
        vjReadLayoutAdaptor: VJReadLayoutAdaptor,
        writer: BufferedWriter)
    {
        // count how many split reads
        val numSplitReads5Bases = overlay.reads.map { o -> vjReadLayoutAdaptor.toReadCandidate(o) }
            .count { o -> Math.max(o.leftSoftClip, o.rightSoftClip) >= 5 }

        val numSplitReads10Bases = overlay.reads.map { o -> vjReadLayoutAdaptor.toReadCandidate(o) }
            .count { o -> Math.max(o.leftSoftClip, o.rightSoftClip) >= 10 }

        val anchorRange = vjReadLayoutAdaptor.getAnchorRange(geneType, overlay)!!
        var sequence = overlay.consensusSequence()
        var support = overlay.highQualSupportString()

        // make sure aligned to codon
        val codonAlignedSeq = sequence.drop(anchorRange.first % 3)
        val aa = insertDashes(Codons.aminoAcidFromBases(codonAlignedSeq), IntRange(anchorRange.first / 3, (anchorRange.last + 1) / 3 - 1))

        // insert - into the sequence and support
        sequence = insertDashes(sequence, anchorRange)
        support = insertDashes(support, anchorRange)

        writer.write("${overlay.id} type: ${geneType}, read count: ${overlay.reads.size}, split(5) read count: ${numSplitReads5Bases}, ")
        writer.write("split(10) read count: ${numSplitReads10Bases}, ")
        writer.write("AA: ${aa}, ")
        writer.write("aligned: ${overlay.alignedPosition}\n")
        writer.write("    ${sequence}\n")
        writer.write("    ${support}\n")

        for (r in overlay.reads)
        {
            val read: VJReadCandidate = vjReadLayoutAdaptor.toReadCandidate(r)
            val readPadding = Math.max(overlay.alignedPosition - r.alignedPosition, 0)
            val paddedSeq = " ".repeat(readPadding) + r.sequence
            val paddedQual = " ".repeat(readPadding) + SAMUtils.phredToFastq(r.baseQualities)
            writer.write("    read: ${read.read}, aligned: ${r.alignedPosition}\n")
            writer.write("    ${insertDashes(paddedSeq, anchorRange)}\n")
            writer.write("    ${insertDashes(paddedQual, anchorRange)}\n")
        }
    }

    private fun insertDashes(str: String, anchorRange: IntRange): String
    {
        return Cdr3Utils.insertDashes(str, anchorRange.first, anchorRange.last + 1)
    }
}