package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.cider.genes.IgTcrConstantDiversityRegion
import com.hartwig.hmftools.cider.genes.IgTcrGeneFile
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import org.apache.logging.log4j.LogManager

object CiderGeneDataLoader
{
    private val sLogger = LogManager.getLogger(CiderGeneDataLoader::class.java)

    fun loadConstantDiversityRegions(refGenomeVersion: RefGenomeVersion): List<IgTcrConstantDiversityRegion>
    {
        val igTcrGenes = IgTcrGeneFile.read(refGenomeVersion)

        val igTcrConstantDiversityRegions = ArrayList<IgTcrConstantDiversityRegion>()
        for (geneData in igTcrGenes)
        {
            if (geneData.geneLocation == null)
                continue

            if (!geneData.geneLocation.inPrimaryAssembly)
                continue

            // we have to fix the naming. It is actually constant + diversity
            if (geneData.region != IgTcrRegion.CONSTANT && geneData.region != IgTcrRegion.D_REGION)
                continue

            val igTcrLocus: IgTcrLocus = IgTcrLocus.fromGeneName(geneData.geneName)

            val constantRegionGene = IgTcrConstantDiversityRegion(igTcrLocus, geneData.geneLocation, geneData.geneName)

            if (igTcrConstantDiversityRegions.contains(constantRegionGene))
                continue

            igTcrConstantDiversityRegions.add(constantRegionGene)
            sLogger.trace("added constant / D region gene: {}, {}",
                geneData.geneName, constantRegionGene)
        }

        return igTcrConstantDiversityRegions
    }

    fun loadAnchorTemplates(refGenomeVersion: RefGenomeVersion): List<VJAnchorTemplate>
    {
        val igTcrGenes = IgTcrGeneFile.read(refGenomeVersion)

        val vjAnchorTemplateList: MutableList<VJAnchorTemplate> = ArrayList()

        for (geneData in igTcrGenes)
        {
            if (geneData.region != IgTcrRegion.V_REGION && geneData.region != IgTcrRegion.J_REGION)
            {
                continue
            }

            if (geneData.anchorSequence == null)
            {
                continue
            }

            val geneType = VJGeneType.fromGeneName(geneData.geneName)

            // we only use primary assembly gene location, need to filter out the ALT locations
            val geneLocation = primaryAssemblyLocation(geneData.geneLocation)
            val anchorLocation = primaryAssemblyLocation(geneData.anchorLocation)

            val vjAnchorTemplate = VJAnchorTemplate(
                geneType, geneData.geneName, geneData.allele, geneLocation, geneData.anchorSequence, anchorLocation)

            vjAnchorTemplateList.add(vjAnchorTemplate)
        }

        return vjAnchorTemplateList
    }

    fun primaryAssemblyLocation(genomicLocation: GenomicLocation?) : GenomicLocation?
    {
        if (genomicLocation == null)
            return null
        return if (genomicLocation.inPrimaryAssembly) genomicLocation else null
    }
}