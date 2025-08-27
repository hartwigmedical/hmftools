package com.hartwig.hmftools.cider.alignment

import com.hartwig.hmftools.cider.CiderFormatter
import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.File

object AlignmentMatchTsvWriter
{
    enum class Column
    {
        cdr3Seq,
        cdr3AA,
        matchType,
        gene,
        functionality,
        pIdent,
        alignStart,
        alignEnd,
        alignScore,
        refStrand,
        refStart,
        refEnd,
        refContig,
        querySeq
    }

    enum class MatchType
    {
        V,
        D,
        J,
        full
    }

    private const val FILE_EXTENSION = ".cider.alignment_match.tsv.gz"

    @JvmStatic
    fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun write(basePath: String, sample: String, alignmentAnnotations: Collection<AlignmentAnnotation>)
    {
        val filePath = generateFilename(basePath, sample)

        val csvFormat = CSVFormat.Builder.create()
            .setDelimiter('\t').setRecordSeparator('\n')
            .setHeader(Column::class.java)
            .setNullString("null")
            .build()

        csvFormat.print(FileWriterUtils.createBufferedWriter(filePath)).use { printer: CSVPrinter ->
            for (ann in alignmentAnnotations)
            {
                writeVDJSequence(printer, ann)
            }
        }
    }

    private fun writeVDJSequence(csvPrinter: CSVPrinter, alignmentAnnotation: AlignmentAnnotation)
    {
        val typeMatch = listOf(Triple(MatchType.V, alignmentAnnotation.vGene, alignmentAnnotation.vMatch),
            Triple(MatchType.D, alignmentAnnotation.dGene, alignmentAnnotation.dMatch),
            Triple(MatchType.J, alignmentAnnotation.jGene, alignmentAnnotation.jMatch),
            Triple(MatchType.full, null, alignmentAnnotation.fullMatch))

        for ((type, gene, match) in typeMatch)
        {
            if (match == null)
                continue

            for (c in Column.values())
            {
                when (c)
                {
                    Column.cdr3Seq -> csvPrinter.print(alignmentAnnotation.vdjSequence.cdr3Sequence)
                    Column.cdr3AA -> csvPrinter.print(CiderFormatter.cdr3AminoAcid(alignmentAnnotation.vdjSequence))
                    Column.matchType -> csvPrinter.print(type)
                    Column.gene -> csvPrinter.print(gene?.geneName)
                    Column.functionality -> csvPrinter.print(gene?.functionality?.toCode())
                    Column.pIdent -> csvPrinter.print(match.percentageIdent)
                    Column.alignStart -> csvPrinter.print(match.queryAlignStart)
                    Column.alignEnd -> csvPrinter.print(match.queryAlignEnd)
                    Column.alignScore -> csvPrinter.print(match.alignmentScore)
                    Column.refStrand -> csvPrinter.print(match.strand.asChar())
                    Column.refStart -> csvPrinter.print(match.refStart)
                    Column.refEnd -> csvPrinter.print(match.refEnd)
                    Column.refContig -> csvPrinter.print(match.refContig)
                    Column.querySeq -> csvPrinter.print(match.querySeq)
                }
            }
            csvPrinter.println()
        }
    }
}