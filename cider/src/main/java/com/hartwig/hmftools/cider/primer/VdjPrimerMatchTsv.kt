package com.hartwig.hmftools.cider.primer

import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.common.utils.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.File

object VdjPrimerMatchTsv
{
    private enum class Column
    {
        cdr3Seq,
        cdr3AA,
        primerGene,
        primerName,
        primerSequence,
        vdjFullSeq
    }

    private const val FILE_EXTENSION = ".cider.primer_match.tsv"

    private fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun writePrimerMatches(outputDir: String, sampleId: String, primerMatches: List<VdjPrimerMatch>)
    {
        val filePath = generateFilename(outputDir, sampleId)

        val csvFormat = CSVFormat.Builder.create()
            .setDelimiter('\t').setRecordSeparator('\n')
            .setHeader(Column::class.java)
            .build()

        csvFormat.print(FileWriterUtils.createBufferedWriter(filePath)).use { writer ->
            for (match in primerMatches)
            {
                writePrimerMatch(writer, match)
            }
        }
    }

    private fun writePrimerMatch(csvPrinter: CSVPrinter, match: VdjPrimerMatch)
    {
        for (c in Column.values())
        {
            when (c)
            {
                Column.cdr3Seq -> csvPrinter.print(match.vdj.cdr3Sequence)
                Column.cdr3AA -> csvPrinter.print(CiderFormatter.cdr3AminoAcid(match.vdj))
                Column.primerGene -> csvPrinter.print(match.primer.target)
                Column.primerName -> csvPrinter.print(match.primer.name)
                Column.primerSequence -> csvPrinter.print(match.primer.sequence)
                Column.vdjFullSeq -> csvPrinter.print(match.fullVdjSequence)
            }
        }
        csvPrinter.println()
    }
}