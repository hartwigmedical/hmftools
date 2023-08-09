package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVRecord
import org.apache.logging.log4j.LogManager
import java.io.BufferedReader
import java.io.InputStreamReader

data class IgTcrGene(
    val geneName: String,
    val chromosome: String,
    val strand: Strand,
    val start: Int,
    val end: Int,
    val geneSegmentType: Char)

object IgTcrGeneLoader
{
    enum class Column
    {
        gene, chromosome, strand, start, end, geneSegmentType
    }

    private val sLogger = LogManager.getLogger(IgTcrGeneLoader::class.java)

    fun load(refGenomeVersion: RefGenomeVersion) : List<IgTcrGene>
    {
        val resourcePath: String = if (refGenomeVersion.is37)
        {
            "igtcr_gene.37.tsv"
        }
        else if (refGenomeVersion.is38)
        {
            "igtcr_gene.38.tsv"
        }
        else
        {
            sLogger.error("unknown ref genome version: {}", refGenomeVersion)
            throw IllegalArgumentException("unknown ref genome version: $refGenomeVersion")
        }

        val igTcrGenes = ArrayList<IgTcrGene>()

        val tsvStream = javaClass.classLoader.getResourceAsStream(resourcePath)
        if (tsvStream == null)
        {
            sLogger.error("unable to find resource file: {}", resourcePath)
            throw RuntimeException("unable to find resource file: $resourcePath")
        }
        BufferedReader(InputStreamReader(tsvStream)).use { reader ->
            val format = CSVFormat.Builder.create()
                .setDelimiter('\t')
                .setRecordSeparator('\n')
                .setHeader().setSkipHeaderRecord(true) // use first line header as column names
                .build()
            val records: Iterable<CSVRecord> = format.parse(reader)
            for (record in records)
            {
                val geneName = record[Column.gene]
                val chromosome = refGenomeVersion.versionedChromosome(record[Column.chromosome])
                val posStart = record[Column.start].toInt()
                val posEnd = record[Column.end].toInt()
                val strand: Strand =  Strand.valueOf(record[Column.strand][0])
                val geneSegmentType: Char = record[Column.geneSegmentType][0]

                igTcrGenes.add(IgTcrGene(geneName = geneName, chromosome=chromosome, strand=strand, start = posStart, end = posEnd,
                    geneSegmentType = geneSegmentType))
            }
        }

        return igTcrGenes
    }
}