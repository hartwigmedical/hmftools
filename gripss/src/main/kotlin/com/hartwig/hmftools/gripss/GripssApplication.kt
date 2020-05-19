package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.gripss.link.AlternatePath
import com.hartwig.hmftools.gripss.link.AssemblyLink
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import htsjdk.variant.vcf.VCFFileReader
import org.apache.logging.log4j.LogManager
import java.io.File

fun main(args: Array<String>) {

    val inputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893R_CPCT02010893T.gridss.vcf.gz"
//    val inputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893T.gridss.somatic.vcf"
    val outputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893T.post.vcf"
    val filterConfig = GripssFilterConfig(
            0.03,
            8,
            0.005,
            0.95,
            1000,
            350,
            50,
            6,
            50,
            5,
            32)
    val config = GripssConfig(inputVCF, outputVCF, filterConfig)

    GripssApplication(config).use { x -> x.run() }
}


class GripssApplication(private val config: GripssConfig) : AutoCloseable, Runnable {

    companion object {
        private val logger = LogManager.getLogger(this::class.java)
    }

    private val startTime = System.currentTimeMillis();
    private val fileReader = VCFFileReader(File(config.inputVcf), false)
    private val fileWriter = GripssVCF(config.outputVcf)

    override fun run() {
        logger.info("Reading file: ${config.inputVcf}")

        fileWriter.writeHeader(fileReader.fileHeader)
        val structuralVariants = hardFilterVariants(fileReader)
        val variantStore = VariantStore.invoke(structuralVariants)

        logger.info("Finding assembly links")
        val assemblyLinks: LinkStore = AssemblyLink(structuralVariants)

        logger.info("Finding transitive links")
        val alternatePaths: Collection<AlternatePath> = AlternatePath(assemblyLinks, variantStore, structuralVariants)
        val transitiveLinks: LinkStore  = LinkStore(alternatePaths.flatMap { x -> x.transitiveLinks() })

        logger.info("Finding double stranded break links")
        // TODO: Add DSB Links
        val allLinks = LinkStore(assemblyLinks, transitiveLinks)

        logger.info("Writing file: ${config.outputVcf}")
        for (variant in structuralVariants) {
            fileWriter.writeVariant(variant.context(config.filterConfig, allLinks.localLinkedBy(variant.vcfId), allLinks.localLinkedBy(variant.mateId), false))
        }

    }

    private fun hardFilterVariants(fileReader: VCFFileReader): List<StructuralVariantContext> {
        val unfiltered: MutableSet<String> = mutableSetOf()
        val hardFilter: MutableSet<String> = mutableSetOf()
        val structuralVariants: MutableList<StructuralVariantContext> = mutableListOf()

        for (variantContext in fileReader) {
            val structuralVariant = StructuralVariantContext(variantContext)
            if (hardFilter.contains(structuralVariant.vcfId) || structuralVariant.isHardFilter(config.filterConfig)) {
                structuralVariant.mateId?.let { hardFilter.add(it) }
            } else {
                unfiltered.add(variantContext.id)
                structuralVariants.add(structuralVariant)
            }
        }

        val mateIsValidOrNull = { x: StructuralVariantContext -> x.mateId?.let { unfiltered.contains(it) } != false}
        return structuralVariants.filter { x -> !hardFilter.contains(x.vcfId) && mateIsValidOrNull(x) }
    }

    override fun close() {
        fileReader.close()
        fileWriter.close()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}

