package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.genes.IgTcrLocus
import com.hartwig.hmftools.cider.genes.VJGeneType
import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.File

object CiderLocusStatsWriter
{
    // locus,readsUsed,readsTotal,downSampled(T/F),sequences,passSequences
    private enum class Column
    {
        locus,
        readsUsed,
        readsTotal,
        downSampled,
        sequences,
        passSequences
    }

    private const val FILE_EXTENSION = ".cider.locus_stats.tsv"

    @JvmStatic
    fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun writeLocusStats(basePath: String, sample: String, layoutBuildResults: Map<VJGeneType, VJReadLayoutBuilder.LayoutBuildResult>,
                        vdjAnnotations: List<VdjAnnotation>)
    {
        val filePath = generateFilename(basePath, sample)

        val csvFormat = CSVFormat.Builder.create()
            .setDelimiter('\t').setRecordSeparator('\n')
            .setHeader(Column::class.java)
            .build()

        csvFormat.print(FileWriterUtils.createBufferedWriter(filePath)).use { printer: CSVPrinter ->
            for (locus in IgTcrLocus.entries)
            {
                writeLocusStats(printer, locus, layoutBuildResults, vdjAnnotations)
            }
        }
    }


    // locus,readsUsed,readsTotal,downSampled(T/F),sequences,passSequences
    private fun writeLocusStats(csvPrinter: CSVPrinter, locus: IgTcrLocus,
                                layoutBuildResults: Map<VJGeneType, VJReadLayoutBuilder.LayoutBuildResult>,
                                vdjAnnotations: List<VdjAnnotation>)
    {
        // gather the data needed for this locus
        val locusLayoutBuildResults = layoutBuildResults.filter { (k, _) -> k.locus == locus }.values

        val readsTotal = locusLayoutBuildResults.sumOf { o -> o.numCandidateReads }

        // to calculate correctly number of reads used, we have to create a set of all
        // the ReadKeys to remove duplicates
        val numReadsUsed: Int = locusLayoutBuildResults
            .flatMap{ o -> o.layouts.flatMap { l -> l.reads } } // list of reads
            .map { r -> r.readKey }.toSet().size

        val downSampled = locusLayoutBuildResults.any { o -> o.downSampled }

        val locusVdjAnnotations = vdjAnnotations.filter { vdjAnnotation -> locus == vdjAnnotation.locus }

        val numSequences = locusVdjAnnotations.size
        val numPassSequences = locusVdjAnnotations.count { vdjAnnotation -> vdjAnnotation.passesFilter }

        for (c in Column.entries)
        {
            when (c)
            {
                Column.locus -> csvPrinter.print(locus.prettyPrint())
                Column.readsUsed -> csvPrinter.print(numReadsUsed)
                Column.readsTotal -> csvPrinter.print(readsTotal)
                Column.downSampled -> csvPrinter.print(downSampled)
                Column.sequences -> csvPrinter.print(numSequences)
                Column.passSequences -> csvPrinter.print(numPassSequences)
            }
        }
        csvPrinter.println()
    }
}