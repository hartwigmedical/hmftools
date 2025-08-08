package com.hartwig.hmftools.cider.blastn

import com.hartwig.hmftools.cider.CiderFormatter
import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.File

object BlastnMatchTsvWriter
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
        refContig
    }

    enum class MatchType
    {
        V,
        D,
        J,
        full
    }

    private const val FILE_EXTENSION = ".cider.blastn_match.tsv.gz"

    @JvmStatic
    fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun write(basePath: String, sample: String, blastnAnnotations: Collection<BlastnAnnotation>)
    {
        val filePath = generateFilename(basePath, sample)

        val csvFormat = CSVFormat.Builder.create()
            .setDelimiter('\t').setRecordSeparator('\n')
            .setHeader(Column::class.java)
            .setNullString("null")
            .build()

        csvFormat.print(FileWriterUtils.createBufferedWriter(filePath)).use { printer: CSVPrinter ->
            for (ann in blastnAnnotations)
            {
                writeVDJSequence(printer, ann)
            }
        }
    }

    private fun writeVDJSequence(csvPrinter: CSVPrinter, blastnAnnotation: BlastnAnnotation)
    {
        val typeMatch = listOf(Triple(MatchType.V, blastnAnnotation.vGene, blastnAnnotation.vMatch),
            Triple(MatchType.D, blastnAnnotation.dGene, blastnAnnotation.dMatch),
            Triple(MatchType.J, blastnAnnotation.jGene, blastnAnnotation.jMatch),
            Triple(MatchType.full, null, blastnAnnotation.fullMatch))

        for ((type, gene, match) in typeMatch)
        {
            if (match == null)
                continue

            for (c in Column.values())
            {
                when (c)
                {
                    Column.cdr3Seq -> csvPrinter.print(blastnAnnotation.vdjSequence.cdr3Sequence)
                    Column.cdr3AA -> csvPrinter.print(CiderFormatter.cdr3AminoAcid(blastnAnnotation.vdjSequence))
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
                }
            }
            csvPrinter.println()
        }
    }
}