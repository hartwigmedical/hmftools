package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.annotation.AlignmentStatus
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
        alignmentStatus,
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
        vGeneSupplementary,
        vPIdent,
        vAlignStart,
        vAlignEnd,
        dGene,
        dGeneSupplementary,
        dPIdent,
        dAlignStart,
        dAlignEnd,
        jGene,
        jGeneSupplementary,
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
                Column.alignmentStatus -> csvPrinter.print(vdjAnnotation.alignmentAnnotation?.alignmentStatus ?: AlignmentStatus.SKIPPED_ALIGN)
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
                Column.vGene -> csvPrinter.print(vdjAnnotation.alignmentAnnotation?.vGene?.geneName)
                Column.vGeneSupplementary -> csvPrinter.print(formatGeneSupplementary(
                    vdjAnnotation.alignmentAnnotation?.vGeneSupplementary, vdjAnnotation.alignmentAnnotation?.vGene))
                Column.vPIdent -> csvPrinter.print(vdjAnnotation.alignmentAnnotation?.vAlignment?.percentageIdent)
                Column.vAlignStart -> if (vdjAnnotation.alignmentAnnotation != null)
                    {
                        csvPrinter.print(zeroBaseAlignStart(vdjAnnotation.alignmentAnnotation!!.vAlignment))
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.vAlignEnd -> if (vdjAnnotation.alignmentAnnotation != null)
                    {
                        csvPrinter.print(vdjAnnotation.alignmentAnnotation!!.vAlignment?.queryEnd)
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.dGene -> csvPrinter.print(vdjAnnotation.alignmentAnnotation?.dGene?.geneName)
                Column.dGeneSupplementary -> csvPrinter.print(formatGeneSupplementary(
                    vdjAnnotation.alignmentAnnotation?.dGeneSupplementary, vdjAnnotation.alignmentAnnotation?.dGene))
                Column.dPIdent -> csvPrinter.print(vdjAnnotation.alignmentAnnotation?.dAlignment?.percentageIdent)
                Column.dAlignStart -> if (vdjAnnotation.alignmentAnnotation != null)
                    {
                        csvPrinter.print(zeroBaseAlignStart(vdjAnnotation.alignmentAnnotation!!.dAlignment))
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.dAlignEnd -> if (vdjAnnotation.alignmentAnnotation != null)
                    {
                        csvPrinter.print(vdjAnnotation.alignmentAnnotation!!.dAlignment?.queryEnd)
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.jGene -> csvPrinter.print(vdjAnnotation.alignmentAnnotation?.jGene?.geneName)
                Column.jGeneSupplementary -> csvPrinter.print(formatGeneSupplementary(
                    vdjAnnotation.alignmentAnnotation?.jGeneSupplementary, vdjAnnotation.alignmentAnnotation?.jGene))
                Column.jPIdent -> csvPrinter.print(vdjAnnotation.alignmentAnnotation?.jAlignment?.percentageIdent)
                Column.jAlignStart -> if (vdjAnnotation.alignmentAnnotation != null)
                    {
                        csvPrinter.print(zeroBaseAlignStart(vdjAnnotation.alignmentAnnotation!!.jAlignment))
                    }
                    else
                    {
                        csvPrinter.print(null)
                    }
                Column.jAlignEnd -> if (vdjAnnotation.alignmentAnnotation != null)
                    {
                        csvPrinter.print(vdjAnnotation.alignmentAnnotation!!.jAlignment?.queryEnd)
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

    private fun zeroBaseAlignStart(alignment: Alignment?) : Int?
    {
        return if (alignment == null) null else alignment.queryStart - 1
    }

    private fun formatGeneSupplementary(geneSupplementary: List<IgTcrGene>?, primaryGene: IgTcrGene?): String?
    {
        return geneSupplementary
            ?.map { it.geneName }
            // Since we only output the gene name, don't output duplicate gene names. Can have multiple alleles with the same name.
            ?.filter { it != primaryGene?.geneName }
            ?.toSet()
            ?.sorted()
            ?.joinToString(";")
    }
}