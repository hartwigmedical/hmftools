package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVRecord
import org.apache.logging.log4j.LogManager
import org.eclipse.collections.impl.list.mutable.FastList
import java.io.BufferedReader
import java.io.IOException
import java.io.InputStreamReader

typealias VdjAnchorTemplateColumn = CiderConstants.VjAnchorTemplateTsvColumn

object CiderGeneDataLoader
{
    private val sLogger = LogManager.getLogger(CiderGeneDataLoader::class.java)

    fun loadConstantRegionGenes(refGenomeVersion: RefGenomeVersion): List<IgTcrConstantRegion>
    {
        val igTcrGenes = IgTcrGeneLoader.load(refGenomeVersion)

        val igTcrConstantRegions = FastList<IgTcrConstantRegion>()
        for (geneData in igTcrGenes)
        {
            if (geneData.geneName.length <= 4) continue
            val geneNamePrefix = geneData.geneName.substring(0, 4)

            val igConstantRegionType: IgTcrConstantRegion.Type = try
            {
                IgTcrConstantRegion.Type.valueOf(geneNamePrefix)
            }
            catch (ignored: IllegalArgumentException)
            {
                continue
            }

            val genomeRegionStrand = GenomeRegionStrand(geneData.chromosome, geneData.start, geneData.end, geneData.strand)

            igTcrConstantRegions.add(IgTcrConstantRegion(igConstantRegionType, genomeRegionStrand))
            sLogger.debug("found constant region gene: {}, type: {}, location: {}",
                geneData.geneName, igConstantRegionType, genomeRegionStrand)
        }

        return igTcrConstantRegions
    }

    @Throws(IOException::class)
    fun loadAnchorTemplateTsv(refGenomeVersion: RefGenomeVersion): List<VJAnchorTemplate>
    {
        val resourcePath: String = if (refGenomeVersion.is37)
        {
            "igtcr_anchor.37.tsv"
        }
        else if (refGenomeVersion.is38)
        {
            "igtcr_anchor.38.tsv"
        }
        else
        {
            sLogger.error("unknown ref genome version: {}", refGenomeVersion)
            throw IllegalArgumentException("unknown ref genome version: $refGenomeVersion")
        }
        val vjAnchorTemplateList: MutableList<VJAnchorTemplate> = ArrayList()
        val tsvStream = CiderGeneDataLoader::class.java.classLoader.getResourceAsStream(resourcePath)
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
                val geneName = record[VdjAnchorTemplateColumn.gene]
                val allele = record[VdjAnchorTemplateColumn.allele]
                val posStart = record[VdjAnchorTemplateColumn.posStart].toInt()
                val posEnd = record[VdjAnchorTemplateColumn.posEnd].toInt()
                val anchorSequence = record[VdjAnchorTemplateColumn.anchorSequence]
                var chromosome = record[VdjAnchorTemplateColumn.chr]
                val strandStr = record[VdjAnchorTemplateColumn.strand]
                var strand: Strand? = null
                if (strandStr == "+") strand = Strand.FORWARD else if (strandStr == "-") strand = Strand.REVERSE
                val anchorStart = record[VdjAnchorTemplateColumn.anchorStart].toInt()
                val anchorEnd = record[VdjAnchorTemplateColumn.anchorEnd].toInt()
                val sequence = record[VdjAnchorTemplateColumn.sequence]
                var genomeRegionStrand: GenomeRegionStrand? = null
                var anchorLocation: GenomeRegionStrand? = null
                if (chromosome.isNotEmpty())
                {
                    chromosome = refGenomeVersion.versionedChromosome(chromosome)
                    if (posStart <= 0 || posEnd <= 0)
                    {
                        throw RuntimeException("chromosome exist but pos start or pos end invalid")
                    }
                    if (strand == null)
                    {
                        throw RuntimeException("chromosome exist but strand invalid")
                    }
                    genomeRegionStrand = GenomeRegionStrand(chromosome, posStart, posEnd, strand)
                    if (anchorStart >= 0 && anchorEnd >= 0)
                    {
                        anchorLocation = GenomeRegionStrand(chromosome, anchorStart, anchorEnd, strand)
                    }
                }

                val geneType = if (geneName == "IGKDEL") VJGeneType.IGKDEL else VJGeneType.valueOf(geneName.substring(0, 4))

                val vjAnchorTemplate = VJAnchorTemplate(
                    geneType, geneName, allele, genomeRegionStrand, sequence, anchorSequence, anchorLocation)
                vjAnchorTemplateList.add(vjAnchorTemplate)
            }
        }
        return vjAnchorTemplateList
    }
}