package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache
import com.hartwig.hmftools.common.gene.GeneData
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVRecord
import org.apache.logging.log4j.LogManager
import org.eclipse.collections.impl.list.mutable.FastList
import java.io.BufferedReader
import java.io.IOException
import java.io.InputStreamReader
import java.util.stream.Collectors

object CiderGeneDataLoader
{
    private val sLogger = LogManager.getLogger(CiderGeneDataLoader::class.java)

    fun loadConstantRegionGenes(refGenomeVersion: RefGenomeVersion, ensemblDataDir: String): List<IgTcrConstantRegion>
    {
        val ensemblDataCache = EnsemblDataCache(ensemblDataDir, refGenomeVersion)
        val ensemblLoadOk = ensemblDataCache.load(true)
        if (!ensemblLoadOk)
        {
            sLogger.error("Ensembl data cache load failed")
            throw RuntimeException("Ensembl data cache load failed")
        }
        val igTcrConstantRegions = FastList<IgTcrConstantRegion>()

        // find all the constant region genes
        val chrGeneDataMap: MutableMap<String, List<GeneData>> = ensemblDataCache.chrGeneDataMap
        val geneDataList = chrGeneDataMap.values.stream().flatMap { obj: List<GeneData> -> obj.stream() }.collect(Collectors.toList())

        // find every gene that is constant region
        for (geneData in geneDataList)
        {
            if (geneData.GeneName.length <= 4) continue
            val geneNamePrefix = geneData.GeneName.substring(0, 4)
            var igConstantRegionType: IgTcrConstantRegion.Type
            igConstantRegionType = try
            {
                IgTcrConstantRegion.Type.valueOf(geneNamePrefix)
            } catch (ignored: IllegalArgumentException)
            {
                continue
            }
            val genomeRegionStrand = GenomeRegionStrand(
                geneData.Chromosome, geneData.GeneStart, geneData.GeneEnd,
                if (geneData.forwardStrand()) Strand.FORWARD else Strand.REVERSE
            )
            igTcrConstantRegions.add(IgTcrConstantRegion(igConstantRegionType, genomeRegionStrand))
            sLogger.info(
                "found constant region gene: {}, type: {}, location: {}",
                geneData.GeneName, igConstantRegionType, genomeRegionStrand
            )
        }
        return igTcrConstantRegions
    }

    @Throws(IOException::class)
    fun loadAnchorTemplateTsv(refGenomeVersion: RefGenomeVersion): List<VJAnchorTemplate>
    {
        val resourcePath: String
        if (refGenomeVersion.is37)
        {
            resourcePath = "igtcr_anchor.37.tsv"
        } else if (refGenomeVersion.is38)
        {
            resourcePath = "igtcr_anchor.38.tsv"
        } else
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
                val id = record["id"]
                val name = record["gene"]
                val allele = record["allele"]
                val posStart = record["posStart"].toInt()
                val posEnd = record["posEnd"].toInt()
                val anchorSequence = record["anchorSequence"]
                var chromosome = record["chr"]
                val strandStr = record["strand"]
                var strand: Strand? = null
                if (strandStr == "+") strand = Strand.FORWARD else if (strandStr == "-") strand = Strand.REVERSE
                val anchorStart = record["anchorStart"].toInt()
                val anchorEnd = record["anchorEnd"].toInt()
                val sequence = record["sequence"]
                var genomeRegionStrand: GenomeRegionStrand? = null
                var anchorLocation: GenomeRegionStrand? = null
                if (!chromosome.isEmpty())
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
                    if (anchorStart >= 0 && anchorEnd >= 0) anchorLocation =
                        GenomeRegionStrand(chromosome, anchorStart, anchorEnd, strand)
                }
                val vjAnchorTemplate = VJAnchorTemplate(
                    id, name, allele, genomeRegionStrand, sequence, anchorSequence, anchorLocation
                )
                vjAnchorTemplateList.add(vjAnchorTemplate)
            }
        }
        return vjAnchorTemplateList
    }
}