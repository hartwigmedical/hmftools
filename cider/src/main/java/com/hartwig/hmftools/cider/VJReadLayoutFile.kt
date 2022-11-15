package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.FileWriterUtils.createGzipBufferedWriter
import htsjdk.samtools.SAMUtils
import java.io.BufferedWriter
import java.io.File

object VJReadLayoutFile
{
    private const val FILE_EXTENSION = ".cider.layout.gz"

    private fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun writeLayouts(outputDir: String, sampleId: String, overlayMap: Map<VJGeneType, List<ReadLayout>>, minSoftClipBases: Int)
    {
        val filePath = generateFilename(outputDir, sampleId)

        createGzipBufferedWriter(filePath).use { writer ->
            for ((geneType, overlayList) in overlayMap)
            {
                writeLayouts(writer, geneType, overlayList)
            }
        }
    }

    private fun writeLayouts(writer: BufferedWriter, geneType: VJGeneType, overlayList: List<ReadLayout>)
    {
        val vjReadLayoutAdaptor = VJReadLayoutBuilder(0, 0)

        // sort the overlays by number of reads
        val sortedOverlayList = overlayList.sortedByDescending({ o -> o.reads.size })

        for (layout in sortedOverlayList)
        {
            writeLayout(layout, geneType, vjReadLayoutAdaptor, writer)
        }
    }

    private fun writeLayout(
        layout: ReadLayout,
        geneType: VJGeneType,
        vjReadLayoutAdaptor: VJReadLayoutBuilder,
        writer: BufferedWriter)
    {
        // count how many split reads
        val numSplitReads5Bases = layout.reads.map { o -> vjReadLayoutAdaptor.toReadCandidate(o) }
            .count { o -> Math.max(o.leftSoftClip, o.rightSoftClip) >= 5 }

        val numSplitReads10Bases = layout.reads.map { o -> vjReadLayoutAdaptor.toReadCandidate(o) }
            .count { o -> Math.max(o.leftSoftClip, o.rightSoftClip) >= 10 }

        val anchorRange = vjReadLayoutAdaptor.getAnchorRange(geneType, layout)

        // if there is no anchor within this layout then don't write it
        if (anchorRange == null)
            return

        var sequence = layout.consensusSequence()
        var support = layout.highQualSupportString()

        // make sure aligned to codon
        val codonAlignedSeq = sequence.drop(anchorRange.first % 3)
        val aa = insertDashes(Codons.aminoAcidFromBases(codonAlignedSeq), IntRange(anchorRange.first / 3, (anchorRange.last + 1) / 3 - 1))

        // insert - into the sequence and support
        sequence = insertDashes(sequence, anchorRange)
        support = insertDashes(support, anchorRange)

        writer.write("${layout.id} type: ${geneType}, read count: ${layout.reads.size}, split(5) read count: ${numSplitReads5Bases}, ")
        writer.write("split(10) read count: ${numSplitReads10Bases}, ")
        writer.write("AA: ${aa}\n")
        writer.write("    ${sequence}\n")
        writer.write("    ${support}\n")

        for (r in layout.reads)
        {
            val read: VJReadCandidate = vjReadLayoutAdaptor.toReadCandidate(r)
            val readPadding = Math.max(layout.alignedPosition - r.alignedPosition, 0)
            val paddedSeq = " ".repeat(readPadding) + r.sequence
            val paddedQual = " ".repeat(readPadding) + SAMUtils.phredToFastq(r.baseQualities)
            writer.write("    read: ${read.read}\n")
            writer.write("    ${insertDashes(paddedSeq, anchorRange)}\n")
            writer.write("    ${insertDashes(paddedQual, anchorRange)}\n")
        }
    }

    private fun insertDashes(str: String, anchorRange: IntRange): String
    {
        return CiderUtils.insertDashes(str, anchorRange.first, anchorRange.last + 1)
    }
}