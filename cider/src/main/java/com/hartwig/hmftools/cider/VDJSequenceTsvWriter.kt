package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.blastn.BlastnStatus
import com.hartwig.hmftools.cider.blastn.BlastnUtil
import com.hartwig.hmftools.common.blastn.BlastnMatch
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.File

object VDJSequenceTsvWriter
{
    enum class Column
    {
        cdr3Seq,
        cdr3AA,
        locus,
        filter,
        blastnStatus,
        minHighQualBaseReads,
        assignedReads,
        vAlignedReads,
        jAlignedReads,
        inFrame,
        containsStop,
        vType,
        vAnchorStart,
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
        jAnchorEnd,
        jAnchorSeq,
        jAnchorTemplateSeq,
        jAnchorAA,
        jAnchorTemplateAA,
        jMatchMethod,
        jSimilarityScore,
        jNonSplitReads,
        vGene,
        vPIdent,
        vAlignStart,
        vAlignEnd,
        dGene,
        dPIdent,
        dAlignStart,
        dAlignEnd,
        jGene,
        jPIdent,
        jAlignStart,
        jAlignEnd,
        vPrimerMatches,
        jPrimerMatches,
        layoutId,
        fullSeq,
        support
    }

    private const val FILE_EXTENSION = ".cider.vdj.tsv.gz"

    @JvmStatic
    fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun writeVDJSequences(basePath: String, sample: String, vdjAnnotations: List<VdjAnnotation>)
    {
        val filePath = generateFilename(basePath, sample)

        val csvFormat = CSVFormat.Builder.create()
            .setDelimiter('\t').setRecordSeparator('\n')
            .setHeader(Column::class.java)
            .build()

        csvFormat.print(FileWriterUtils.createGzipBufferedWriter(filePath)).use { printer: CSVPrinter ->
            for (vdjAnn in vdjAnnotations)
            {
                writeVDJSequence(printer, vdjAnn)
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
                Column.locus -> csvPrinter.print(vdjAnnotation.locus.prettyPrint())
                Column.filter -> csvPrinter.print(vdjAnnotation.filters.joinToString(separator = ";"))
                Column.blastnStatus -> csvPrinter.print(vdjAnnotation.blastnAnnotation?.blastnStatus ?: BlastnStatus.SKIPPED_BLASTN)
                Column.minHighQualBaseReads -> csvPrinter.print(vdjAnnotation.cdr3SupportMin)
                Column.assignedReads -> csvPrinter.print(vdj.numReads)
                Column.vAlignedReads -> csvPrinter.print(vdjAnnotation.vAlignedReads)
                Column.jAlignedReads -> csvPrinter.print(vdjAnnotation.jAlignedReads)
                Column.inFrame -> csvPrinter.print(vdj.isInFrame)
                Column.containsStop -> csvPrinter.print(vdj.aminoAcidSequence.contains(Codons.STOP_AMINO_ACID))
                Column.vType -> csvPrinter.print(vdj.vAnchor?.geneType)
                Column.vAnchorStart -> csvPrinter.print(if (vdj.vAnchor != null) vdj.layoutSliceStart else null)
                Column.vAnchorEnd -> csvPrinter.print(if (vdj.vAnchor != null) vdj.layoutSliceStart + vdj.vAnchor.anchorBoundary else null)
                Column.vAnchorSeq -> csvPrinter.print(vdj.vAnchorSequence)
                Column.vAnchorTemplateSeq -> csvPrinter.print(vdj.vAnchor?.templateAnchorSeq)
                Column.vAnchorAA -> csvPrinter.print(CiderFormatter.vAnchorAA(vdj))
                Column.vAnchorTemplateAA -> csvPrinter.print(if (vdj.vAnchor != null) CiderFormatter.aminoAcidFromBases(vdj.vAnchor.templateAnchorSeq) else null)
                Column.vMatchMethod -> csvPrinter.print(vdj.vAnchor?.matchMethod)
                Column.vSimilarityScore -> csvPrinter.print(vdjAnnotation.vSimilarityScore)
                Column.vNonSplitReads -> csvPrinter.print(vdjAnnotation.vNonSplitReads)
                Column.jType -> csvPrinter.print(vdj.jAnchor?.geneType)
                Column.jAnchorStart -> csvPrinter.print(if (vdj.jAnchor != null) vdj.layoutSliceStart + vdj.jAnchor.anchorBoundary else null)
                Column.jAnchorEnd -> csvPrinter.print(if (vdj.jAnchor != null) vdj.layoutSliceEnd else null)
                Column.jAnchorSeq -> csvPrinter.print(vdj.jAnchorSequence)
                Column.jAnchorTemplateSeq -> csvPrinter.print(vdj.jAnchor?.templateAnchorSeq)
                Column.jAnchorAA -> csvPrinter.print(CiderFormatter.jAnchorAA(vdj))
                Column.jAnchorTemplateAA -> csvPrinter.print(if (vdj.jAnchor != null) CiderFormatter.aminoAcidFromBases(vdj.jAnchor.templateAnchorSeq) else null)
                Column.jMatchMethod -> csvPrinter.print(vdj.jAnchor?.matchMethod)
                Column.jSimilarityScore -> csvPrinter.print(vdjAnnotation.jSimilarityScore)
                Column.jNonSplitReads -> csvPrinter.print(vdjAnnotation.jNonSplitReads)
                Column.vGene -> csvPrinter.print(vdjAnnotation.blastnAnnotation?.vGene?.geneName)
                Column.vPIdent -> csvPrinter.print(vdjAnnotation.blastnAnnotation?.vMatch?.percentageIdent)
                Column.vAlignStart -> if (vdjAnnotation.blastnAnnotation != null)
                    {
                        csvPrinter.print(zeroBaseAlignStart(vdjAnnotation.blastnAnnotation!!.vMatch))
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.vAlignEnd -> if (vdjAnnotation.blastnAnnotation != null)
                    {
                        csvPrinter.print(vdjAnnotation.blastnAnnotation!!.vMatch?.queryAlignEnd)
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.dGene -> csvPrinter.print(vdjAnnotation.blastnAnnotation?.dGene?.geneName)
                Column.dPIdent -> csvPrinter.print(vdjAnnotation.blastnAnnotation?.dMatch?.percentageIdent)
                Column.dAlignStart -> if (vdjAnnotation.blastnAnnotation != null)
                    {
                        csvPrinter.print(zeroBaseAlignStart(vdjAnnotation.blastnAnnotation!!.dMatch))
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.dAlignEnd -> if (vdjAnnotation.blastnAnnotation != null)
                    {
                        csvPrinter.print(vdjAnnotation.blastnAnnotation!!.dMatch?.queryAlignEnd)
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.jGene -> csvPrinter.print(vdjAnnotation.blastnAnnotation?.jGene?.geneName)
                Column.jPIdent -> csvPrinter.print(vdjAnnotation.blastnAnnotation?.jMatch?.percentageIdent)
                Column.jAlignStart -> if (vdjAnnotation.blastnAnnotation != null)
                    {
                        csvPrinter.print(zeroBaseAlignStart(vdjAnnotation.blastnAnnotation!!.jMatch))
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.jAlignEnd -> if (vdjAnnotation.blastnAnnotation != null)
                    {
                        csvPrinter.print(vdjAnnotation.blastnAnnotation!!.jMatch?.queryAlignEnd)
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.vPrimerMatches -> csvPrinter.print(vdjAnnotation.vPrimerMatchCount)
                Column.jPrimerMatches -> csvPrinter.print(vdjAnnotation.jPrimerMatchCount)
                Column.layoutId -> csvPrinter.print(vdj.layout.id)
                Column.fullSeq -> csvPrinter.print(vdj.layout.consensusSequenceString())
                Column.support -> csvPrinter.print(CiderUtils.countsToString(vdj.layout.highQualSupportCounts()))
            }
        }
        csvPrinter.println()
    }

    private fun zeroBaseAlignStart(blastnMatch: BlastnUtil.BwaMemMatch?) : Int?
    {
        return if (blastnMatch == null) null else blastnMatch.queryAlignStart - 1
    }
}