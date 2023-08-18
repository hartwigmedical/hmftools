package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.genes.IgTcrConstantRegion
import com.hartwig.hmftools.cider.genes.IgTcrGeneFile
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import org.apache.logging.log4j.LogManager

object CiderGeneDataLoader
{
    private val sLogger = LogManager.getLogger(CiderGeneDataLoader::class.java)

    fun loadConstantRegionGenes(refGenomeVersion: RefGenomeVersion): List<IgTcrConstantRegion>
    {
        val igTcrGenes = IgTcrGeneFile.read(refGenomeVersion)

        val igTcrConstantRegions = ArrayList<IgTcrConstantRegion>()
        for (geneData in igTcrGenes)
        {
            if (geneData.geneName.length <= 4) continue

            if (geneData.geneLocation == null)
                continue

            if (!geneData.geneLocation.isPrimaryAssembly)
                continue

            if (geneData.region != IgTcrRegion.CONSTANT)
                continue

            val igTcrLocus: IgTcrLocus = try
            {
                IgTcrLocus.valueOf(geneData.geneName.substring(0, 3))
            }
            catch (ignored: IllegalArgumentException)
            {
                continue
            }

            val constantRegionGene = IgTcrConstantRegion(igTcrLocus, geneData.geneLocation)

            if (igTcrConstantRegions.contains(constantRegionGene))
                continue

            igTcrConstantRegions.add(constantRegionGene)
            sLogger.debug("added constant region gene: {}, {}",
                geneData.geneName, constantRegionGene)
        }

        return igTcrConstantRegions
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

            val vjAnchorTemplate = VJAnchorTemplate(
                geneType, geneData.geneName, geneData.allele, geneData.geneLocation, geneData.anchorSequence, geneData.anchorLocation)

            vjAnchorTemplateList.add(vjAnchorTemplate)
        }

        return vjAnchorTemplateList
    }
}