package com.hartwig.hmftools.cider.primer

import com.hartwig.hmftools.cider.VDJSequence
import com.hartwig.hmftools.cider.VJReadCandidate
import com.hartwig.hmftools.cider.VJReadLayoutAdaptor
import com.hartwig.hmftools.common.utils.FileWriterUtils.createGzipBufferedWriter
import htsjdk.samtools.SAMUtils
import java.io.BufferedWriter
import java.io.File

object FullSequenceLayoutFile
{
    private const val FILE_EXTENSION = ".cider.full_seq.gz"

    private fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun writeLayouts(outputDir: String, sampleId: String, vdjs: List<VDJSequence>, vjReadLayoutAdaptor: VJReadLayoutAdaptor)
    {
        val filePath = generateFilename(outputDir, sampleId)

        createGzipBufferedWriter(filePath).use { writer ->
            for (vdj in vdjs)
            {
                writeLayout(writer, vdj, vjReadLayoutAdaptor)
            }
        }
    }

    private fun writeLayout(
        writer: BufferedWriter,
        vdj: VDJSequence,
        vjReadLayoutAdaptor: VJReadLayoutAdaptor)
    {
        val layout = vjReadLayoutAdaptor.buildFullLayout(vdj.layout)
        var sequence = layout.consensusSequence()
        var support = layout.highQualSupportString()

        writer.write("${layout.id} read count: ${layout.reads.size}, vdj: ${vdj.aminoAcidSequenceFormatted}")
        writer.write(", aligned: ${layout.alignedPosition}\n")
        writer.write("    ${sequence}\n")
        writer.write("    ${support}\n")

        for (r in layout.reads)
        {
            val read: VJReadCandidate = vjReadLayoutAdaptor.toReadCandidate(r)
            val readPadding = Math.max(layout.alignedPosition - r.alignedPosition, 0)
            val paddedSeq = " ".repeat(readPadding) + r.sequence
            val paddedQual = " ".repeat(readPadding) + SAMUtils.phredToFastq(r.baseQualities)
            writer.write("    read: ${read.read}, aligned: ${r.alignedPosition}\n")
            writer.write("    ${paddedSeq}\n")
            writer.write("    ${paddedQual}\n")
        }
    }
}