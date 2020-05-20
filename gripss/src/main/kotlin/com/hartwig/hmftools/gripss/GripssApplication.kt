package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.bedpe.dedup.DedupPair
import com.hartwig.hmftools.bedpe.dedup.DedupSingle
import com.hartwig.hmftools.gripss.link.AlternatePath
import com.hartwig.hmftools.gripss.link.AssemblyLink
import com.hartwig.hmftools.gripss.link.DsbLink
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore
import htsjdk.variant.vcf.VCFFileReader
import org.apache.logging.log4j.LogManager
import java.io.File

fun main(args: Array<String>) {

    val inputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893R_CPCT02010893T.gridss.vcf.gz"
    val outputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893T.post.vcf"
    val filterConfig = GripssFilterConfig.default()
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

        logger.info("Initial soft filters")
        val softFiltersInitial = SoftFilterStore(config.filterConfig, structuralVariants)

        logger.info("Finding assembly links")
        val assemblyLinks: LinkStore = AssemblyLink(structuralVariants)

        logger.info("Finding transitive links")
        val alternatePaths: Collection<AlternatePath> = AlternatePath(assemblyLinks, variantStore, structuralVariants)
        val alternatePathsStringsByVcfId = alternatePaths.associate { x -> Pair(x.vcfId, x.pathString())}
        val transitiveLinks = LinkStore(alternatePaths.flatMap { x -> x.transitiveLinks() })

        logger.info("Paired break end de-duplication")
        val dedupPair = DedupPair(softFiltersInitial, alternatePaths)
        val softFiltersAfterPairedDedup = softFiltersInitial.update(dedupPair.duplicates, dedupPair.rescue)

        logger.info("Single break end de-duplication")
        val dedupSingle = DedupSingle(variantStore, softFiltersAfterPairedDedup, structuralVariants)
        val softFiltersAfterSingleDedup = softFiltersAfterPairedDedup.update(dedupSingle.duplicates, setOf())

        logger.info("Finding double stranded break links")
        val dsbLinks = DsbLink(variantStore, assemblyLinks, softFiltersAfterSingleDedup.duplicates(), structuralVariants)
        val combinedLinks = LinkStore(assemblyLinks, transitiveLinks, dsbLinks)

        logger.info("Writing file: ${config.outputVcf}")
        for (variant in structuralVariants) {

            val localLinkedBy = combinedLinks[variant.vcfId]
            val remoteLinkedBy = combinedLinks[variant.mateId]
            val filters = softFiltersAfterSingleDedup[variant.vcfId]
            val altPath = alternatePathsStringsByVcfId[variant.vcfId]


            fileWriter.writeVariant(variant.context(localLinkedBy, remoteLinkedBy, altPath, filters))
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

        val mateIsValidOrNull = { x: StructuralVariantContext -> x.mateId?.let { unfiltered.contains(it) } != false }
        return structuralVariants.filter { x -> !hardFilter.contains(x.vcfId) && mateIsValidOrNull(x) }
    }

    override fun close() {
        fileReader.close()
        fileWriter.close()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}

