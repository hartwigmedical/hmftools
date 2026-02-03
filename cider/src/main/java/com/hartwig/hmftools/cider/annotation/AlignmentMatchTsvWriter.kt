package com.hartwig.hmftools.cider.annotation

import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.File

object AlignmentMatchTsvWriter
{
    enum class Column
    {
        cdr3Seq,
        fullSeq,
        matchType,
        gene,
        functionality,
        layoutAlignStart,
        layoutAlignEnd,
        alignScore,
        strand,
        refStart,
        refEnd,
        refContig,
        refContigLength,
        cigar,
        editDistance,
        querySeqStart,
        querySeqEnd,
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
        val typeMatch = listOf(
            Triple(MatchType.V, alignmentAnnotation.vGene?.gene, alignmentAnnotation.vGene?.alignment),
            Triple(MatchType.D, alignmentAnnotation.dGene?.gene, alignmentAnnotation.dGene?.alignment),
            Triple(MatchType.J, alignmentAnnotation.jGene?.gene, alignmentAnnotation.jGene?.alignment),
            Triple(MatchType.full, null, alignmentAnnotation.fullAlignment))

        for ((type, gene, alignment) in typeMatch)
        {
            if (alignment == null)
                continue

            val alignmentOffset = alignmentAnnotation.alignmentQueryRange.start

            for (c in Column.entries)
            {
                when (c)
                {
                    Column.cdr3Seq -> csvPrinter.print(alignmentAnnotation.vdjSequence.cdr3Sequence)
                    Column.fullSeq -> csvPrinter.print(alignmentAnnotation.vdjSequence.layout.consensusSequenceString())
                    Column.matchType -> csvPrinter.print(type)
                    Column.gene -> csvPrinter.print(gene?.geneAllele)
                    Column.functionality -> csvPrinter.print(gene?.functionality?.toCode())
                    Column.layoutAlignStart -> csvPrinter.print(alignmentOffset + alignment.queryStart)
                    Column.layoutAlignEnd -> csvPrinter.print(alignmentOffset + alignment.queryEnd)
                    Column.alignScore -> csvPrinter.print(alignment.alignmentScore)
                    Column.strand -> csvPrinter.print(alignment.strand.asChar())
                    Column.refStart -> csvPrinter.print(alignment.refStart)
                    Column.refEnd -> csvPrinter.print(alignment.refEnd)
                    Column.refContig -> csvPrinter.print(alignment.refContig)
                    Column.refContigLength -> csvPrinter.print(alignment.refContigLength)
                    Column.cigar -> csvPrinter.print(alignment.cigar.joinToString("") { it.toString() })
                    Column.editDistance -> csvPrinter.print(alignment.editDistance)
                    Column.querySeqStart -> csvPrinter.print(alignmentAnnotation.alignmentQueryRange.start)
                    Column.querySeqEnd -> csvPrinter.print(alignmentAnnotation.alignmentQueryRange.endInclusive + 1)
                    Column.querySeq -> csvPrinter.print(alignment.querySeq)
                }
            }
            csvPrinter.println()
        }
    }
}