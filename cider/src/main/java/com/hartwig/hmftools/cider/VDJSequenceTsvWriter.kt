package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.File

object VDJSequenceTsvWriter
{
    private enum class Column
    {
        cdr3Seq,
        cdr3AA,
        filter,
        minHighQualBaseReads,
        assignedReads,
        vAlignedReads,
        jAlignedReads,
        inFrame,
        containsStop,
        vType,
        vAnchorEnd,
        vAnchorSeq,
        vAnchorTemplateSeq,
        vAnchorAA,
        vAnchorTemplateAA,
        vMatchMethod,
        vSimilarityScore,
        vNonSplitReads,
        jType,
        jAnchorStart,
        jAnchorSeq,
        jAnchorTemplateSeq,
        jAnchorAA,
        jAnchorTemplateAA,
        jMatchMethod,
        jSimilarityScore,
        jNonSplitReads,
        vPrimerMatches,
        jPrimerMatches,
        layoutId,
        vdjSeq,
        support
    }

    private const val FILE_EXTENSION = ".cider.vdj_seq.tsv"

    @JvmStatic
    fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun writeVDJSequences(
        basePath: String, sample: String, vdjAnnotations: List<VdjAnnotation>, reportPartialSeq: Boolean)
    {
        val filePath = generateFilename(basePath, sample)

        val csvFormat = CSVFormat.Builder.create()
            .setDelimiter('\t').setRecordSeparator('\n')
            .setHeader(Column::class.java)
            .build()

        csvFormat.print(FileWriterUtils.createBufferedWriter(filePath)).use { printer: CSVPrinter ->
            for (vdjAnn in vdjAnnotations)
            {
                if (reportPartialSeq || vdjAnn.vdj.isFullyRearranged)
                {
                    writeVDJSequence(printer, vdjAnn)
                }
            }
        }
    }

    private fun writeVDJSequence(csvPrinter: CSVPrinter, vdjAnnotation: VdjAnnotation)
    {
        val vdj = vdjAnnotation.vdj

        for (c in Column.values())
        {
            when (c)
            {
                Column.cdr3Seq -> csvPrinter.print(vdj.cdr3Sequence)
                Column.cdr3AA -> csvPrinter.print(CiderFormatter.cdr3AminoAcid(vdj))
                Column.filter -> csvPrinter.print(vdjAnnotation.filters.joinToString(separator = ";"))
                Column.minHighQualBaseReads -> csvPrinter.print(vdj.cdr3SupportMin)
                Column.assignedReads -> csvPrinter.print(vdj.numReads)
                Column.vAlignedReads -> csvPrinter.print(vdjAnnotation.vAlignedReads)
                Column.jAlignedReads -> csvPrinter.print(vdjAnnotation.jAlignedReads)
                Column.inFrame -> csvPrinter.print(vdj.isInFrame)
                Column.containsStop -> csvPrinter.print(vdj.aminoAcidSequence.contains(Codons.STOP_AMINO_ACID))
                Column.vType -> csvPrinter.print(vdj.vAnchor?.geneType)
                Column.vAnchorEnd -> csvPrinter.print(vdj.vAnchor?.anchorBoundary)
                Column.vAnchorSeq -> csvPrinter.print(vdj.vAnchorSequence)
                Column.vAnchorTemplateSeq -> csvPrinter.print(vdj.vAnchor?.templateAnchorSeq)
                Column.vAnchorAA -> csvPrinter.print(CiderFormatter.vAnchorAA(vdj))
                Column.vAnchorTemplateAA -> csvPrinter.print(if (vdj.vAnchor != null) CiderFormatter.aminoAcidFromBases(vdj.vAnchor.templateAnchorSeq) else null)
                Column.vMatchMethod -> csvPrinter.print(vdj.vAnchor?.matchMethod)
                Column.vSimilarityScore -> csvPrinter.print(vdjAnnotation.vSimilarityScore)
                Column.vNonSplitReads -> csvPrinter.print(vdjAnnotation.vNonSplitReads)
                Column.jType -> csvPrinter.print(vdj.jAnchor?.geneType)
                Column.jAnchorStart -> csvPrinter.print(vdj.jAnchor?.anchorBoundary)
                Column.jAnchorSeq -> csvPrinter.print(vdj.jAnchorSequence)
                Column.jAnchorTemplateSeq -> csvPrinter.print(vdj.jAnchor?.templateAnchorSeq)
                Column.jAnchorAA -> csvPrinter.print(CiderFormatter.jAnchorAA(vdj))
                Column.jAnchorTemplateAA -> csvPrinter.print(if (vdj.jAnchor != null) CiderFormatter.aminoAcidFromBases(vdj.jAnchor.templateAnchorSeq) else null)
                Column.jMatchMethod -> csvPrinter.print(vdj.jAnchor?.matchMethod)
                Column.jSimilarityScore -> csvPrinter.print(vdjAnnotation.jSimilarityScore)
                Column.jNonSplitReads -> csvPrinter.print(vdjAnnotation.jNonSplitReads)
                Column.vPrimerMatches -> csvPrinter.print(vdjAnnotation.vPrimerMatchCount)
                Column.jPrimerMatches -> csvPrinter.print(vdjAnnotation.jPrimerMatchCount)
                Column.layoutId -> csvPrinter.print(vdj.layout.id)
                Column.vdjSeq -> csvPrinter.print(vdj.sequence)
                Column.support -> csvPrinter.print(CiderUtils.countsToString(vdj.supportCounts))
            }
        }
        csvPrinter.println()
    }
}