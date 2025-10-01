package com.hartwig.hmftools.cider.annotation

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
        val typeMatch = listOf(Triple(MatchType.V, alignmentAnnotation.vGene, alignmentAnnotation.vAlignment),
            Triple(MatchType.D, alignmentAnnotation.dGene, alignmentAnnotation.dAlignment),
            Triple(MatchType.J, alignmentAnnotation.jGene, alignmentAnnotation.jAlignment),
            Triple(MatchType.full, null, alignmentAnnotation.fullAlignment))

        for ((type, gene, alignment) in typeMatch)
        {
            if (alignment == null)
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
                    Column.pIdent -> csvPrinter.print(alignment.percentageIdent)
                    Column.alignStart -> csvPrinter.print(alignment.queryAlignStart)
                    Column.alignEnd -> csvPrinter.print(alignment.queryAlignEnd)
                    Column.alignScore -> csvPrinter.print(alignment.alignmentScore)
                    Column.refStrand -> csvPrinter.print(alignment.refStrand.asChar())
                    Column.refStart -> csvPrinter.print(alignment.refStart)
                    Column.refEnd -> csvPrinter.print(alignment.refEnd)
                    Column.refContig -> csvPrinter.print(alignment.refContig)
                    Column.querySeq -> csvPrinter.print(alignment.querySeq)
                }
            }
            csvPrinter.println()
        }
    }
}