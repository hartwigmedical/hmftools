package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.genes.VJGeneType
import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.file.FileWriterUtils.createGzipBufferedWriter
import htsjdk.samtools.SAMUtils
import java.io.BufferedWriter
import java.io.File

object VJReadLayoutFile
{
    private const val FILE_EXTENSION = ".cider.layout.gz"

    const val ANSI_RESET = "\u001B[0m"
    const val ANSI_RED = "\u001B[31m"
    const val ANSI_GREEN = "\u001B[32m"
    const val ANSI_YELLOW = "\u001B[33m"
    const val ANSI_PURPLE = "\u001B[35m"

    private fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun writeLayouts(outputDir: String, sampleId: String, overlayMap: Map<VJGeneType, List<ReadLayout>>)
    {
        val filePath = generateFilename(outputDir, sampleId)

        createGzipBufferedWriter(filePath).use { writer ->
            for ((geneType, overlayList) in overlayMap.entries.sortedBy { it.key })
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

        val anchorRange = vjReadLayoutAdaptor.getExtrapolatedAnchorRange(geneType.vj, layout)

        var sequence = layout.consensusSequenceString()
        var support = layout.highQualSupportString()

        // make sure aligned to codon
        val codonAlignedSeq = sequence.drop(Math.floorMod(anchorRange.first, 3))
        val aa = insertDashes(Codons.aminoAcidFromBases(codonAlignedSeq), IntRange(anchorRange.first / 3, (anchorRange.last + 1) / 3 - 1))

        // insert - into the sequence and support
        sequence = coloriseSequence(insertDashes(sequence, anchorRange))
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
            val paddedSeq = " ".repeat(readPadding) + r.sequenceString
            val paddedQual = " ".repeat(readPadding) + SAMUtils.phredToFastq(r.baseQualities)
            writer.write("    read: ${read.read}\n")
            writer.write("    ${coloriseSequence(insertDashes(paddedSeq, anchorRange))}\n")
            writer.write("    ${insertDashes(paddedQual, anchorRange)}\n")
        }
    }

    private fun insertDashes(str: String, anchorRange: IntRange): String
    {
        return CiderUtils.insertDashes(str, anchorRange.first, anchorRange.last + 1)
    }

    private fun coloriseSequence(str: String): String
    {
        // use IGV colours
        // adenine in green, cytosine in blue, guanine in yellow, and thymine in red (A, C, G, T).
        val stringBuilder = StringBuilder()

        for (b in str)
        {
            when (b)
            {
                'A' -> stringBuilder.append(ANSI_GREEN)
                'C' -> stringBuilder.append(ANSI_PURPLE)
                'G' -> stringBuilder.append(ANSI_YELLOW)
                'T' -> stringBuilder.append(ANSI_RED)
            }

            stringBuilder.append(b)
        }

        stringBuilder.append(ANSI_RESET)
        return stringBuilder.toString()
    }
}