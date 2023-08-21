package com.hartwig.hmftools.cider.genes

import com.hartwig.hmftools.cider.IgTcrFunctionality
import com.hartwig.hmftools.cider.IgTcrGene
import com.hartwig.hmftools.cider.IgTcrRegion
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVRecord
import org.apache.logging.log4j.LogManager
import java.io.BufferedReader
import java.io.InputStreamReader

object IgTcrGeneFile
{
    enum class Column
    {
        gene,
        allele,
        region,
        functionality,
        primaryAssembly,
        assemblyName,
        chromosome,
        posStart,
        posEnd,
        strand,
        anchorStart,
        anchorEnd,
        anchorSequence,
        anchorAA
    }

    private val sLogger = LogManager.getLogger(IgTcrGeneFile::class.java)

    fun read(refGenomeVersion: RefGenomeVersion) : List<IgTcrGene>
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
                val allele = record[Column.allele]
                val region = IgTcrRegion.valueOf(record[Column.region])
                val functionality = IgTcrFunctionality.fromCode(record[Column.functionality])
                val isPrimaryAssembly = record[Column.primaryAssembly].toBoolean()
                val assemblyName = if (isPrimaryAssembly) null else record[Column.assemblyName].intern()
                var anchorSequence: String? = record[Column.anchorSequence]
                if (anchorSequence!!.isEmpty())
                    anchorSequence = null
                var chromosome = record[Column.chromosome].intern()
                var genomicLocation: GenomicLocation? = null
                var anchorLocation: GenomicLocation? = null
                if (chromosome.isNotEmpty())
                {
                    chromosome = refGenomeVersion.versionedChromosome(chromosome)

                    val posStart = record[Column.posStart].toInt()
                    val posEnd = record[Column.posEnd].toInt()

                    val strandStr = record[Column.strand]

                    if (strandStr.isEmpty())
                    {
                        throw RuntimeException("chromosome exist but strand is missing")
                    }

                    val strand: Strand = Strand.fromChar(strandStr[0])

                    if (posStart <= 0 || posEnd <= 0)
                    {
                        throw RuntimeException("chromosome exist but pos start or pos end invalid")
                    }
                    genomicLocation = GenomicLocation(chromosome, posStart, posEnd, strand, assemblyName)

                    val anchorStart = record[Column.anchorStart]
                    val anchorEnd = record[Column.anchorEnd]

                    if (anchorStart.isNotEmpty() && anchorEnd.isNotEmpty())
                    {
                        anchorLocation = GenomicLocation(chromosome, anchorStart.toInt(), anchorEnd.toInt(), strand, assemblyName)
                    }
                }
                
                igTcrGenes.add(
                    IgTcrGene(geneName = geneName, allele = allele, region = region, functionality = functionality,
                    geneLocation = genomicLocation, anchorSequence = anchorSequence, anchorLocation = anchorLocation)
                )
            }
        }

        return igTcrGenes
    }

    fun write(tsvPath: String, geneList: List<IgTcrGene>)
    {
        val csvFormat = CSVFormat.Builder.create()
            .setDelimiter('\t').setRecordSeparator('\n')
            .setHeader(Column::class.java)
            .build()

        csvFormat.print(FileWriterUtils.createBufferedWriter(tsvPath)).use { csvPrinter ->
            for (gene: IgTcrGene in geneList)
            {
                for (col in Column.values())
                {
                    when (col)
                    {
                        Column.gene -> csvPrinter.print(gene.geneName)
                        Column.allele -> csvPrinter.print(gene.allele)
                        Column.region -> csvPrinter.print(gene.region)
                        Column.functionality -> csvPrinter.print(gene.functionality.toCode())
                        Column.primaryAssembly -> csvPrinter.print(gene.geneLocation?.inPrimaryAssembly)
                        Column.assemblyName -> csvPrinter.print(gene.geneLocation?.altAssemblyName)
                        Column.chromosome -> csvPrinter.print(gene.geneLocation?.chromosome)
                        Column.posStart -> csvPrinter.print(gene.geneLocation?.posStart)
                        Column.posEnd -> csvPrinter.print(gene.geneLocation?.posEnd)
                        Column.strand -> csvPrinter.print(gene.geneLocation?.strand?.asChar())
                        Column.anchorStart -> csvPrinter.print(gene.anchorLocation?.posStart)
                        Column.anchorEnd -> csvPrinter.print(gene.anchorLocation?.posEnd)
                        Column.anchorSequence -> csvPrinter.print(gene.anchorSequence)
                        Column.anchorAA -> csvPrinter.print(if (gene.anchorSequence != null) Codons.aminoAcidFromBases(gene.anchorSequence) else null)
                        //Column.sequence -> csvPrinter.print(gene.sequence)
                    }
                }
                csvPrinter.println()
            }
        }
    }
}