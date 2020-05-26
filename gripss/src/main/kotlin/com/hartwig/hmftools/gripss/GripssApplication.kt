package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.bedpe.Breakend
import com.hartwig.hmftools.bedpe.Breakpoint
import com.hartwig.hmftools.gripss.dedup.DedupPair
import com.hartwig.hmftools.gripss.dedup.DedupSingle
import com.hartwig.hmftools.gripss.link.AlternatePath
import com.hartwig.hmftools.gripss.link.AssemblyLink
import com.hartwig.hmftools.gripss.link.DsbLink
import com.hartwig.hmftools.gripss.link.LinkRescue
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.LocationStore
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.variant.vcf.VCFFileReader
import org.apache.logging.log4j.LogManager
import java.io.File

fun main(args: Array<String>) {

    val singlePonFile = "/Users/jon/hmf/resources/gridss_pon_single_breakend.bed"
    val pairedPonFile = "/Users/jon/hmf/resources/gridss_pon_breakpoint.bedpe"
    val pairedHotspotFile = "/Users/jon/hmf/resources/gridss_hotspot_breakpoint.bedpe"
    val inputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893R_CPCT02010893T.gridss.vcf.gz"
//    val inputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893R_CPCT02010893T.gridss.chr1.vcf"
    val outputVCF = "/Users/jon/hmf/analysis/gridss/CPCT02010893T.post.vcf"
    val refGenome = "/Users/jon/hmf/resources/Homo_sapiens.GRCh37.GATK.illumina.fasta"
    val filterConfig = GripssFilterConfig.default()
    val config = GripssConfig(inputVCF, outputVCF, singlePonFile, pairedPonFile, pairedHotspotFile, refGenome, filterConfig)

    GripssApplication(config).use { x -> x.run() }
}


class GripssApplication(private val config: GripssConfig) : AutoCloseable, Runnable {
    companion object {
        private val logger = LogManager.getLogger(this::class.java)
        const val PON_ADDITIONAL_DISTANCE = 0
    }

    private val startTime = System.currentTimeMillis();
    private val fileReader = VCFFileReader(File(config.inputVcf), false)
    private val fileWriter = GripssVCF(config.outputVcf)
    private val refGenome = IndexedFastaSequenceFile(File(config.refGenome))

    override fun run() {
        logger.info("Reading VCF file: ${config.inputVcf}")
        val contigComparator = ContigComparator(fileReader.fileHeader.sequenceDictionary)
        val variantStore = VariantStore(hardFilterVariants(fileReader))

        logger.info("Reading hotspot file: ${config.pairedHotspotFile}")
        val hotspotStore = LocationStore(listOf(), Breakpoint.fromBedpeFile(config.pairedHotspotFile, contigComparator))

        logger.info("Reading PON files: ${config.singlePonFile} ${config.pairedPonFile}")
        val breakends = Breakend.fromBedFile(config.singlePonFile)
        val breakpoints = Breakpoint.fromBedpeFile(config.pairedPonFile, contigComparator)
        val ponStore = LocationStore(breakends, breakpoints, PON_ADDITIONAL_DISTANCE)
//        val ponStore = LocationStore(listOf(), listOf(), PON_ADDITIONAL_DISTANCE)



        logger.info("Identifying hotspot variants")
        val hotspots = hotspotStore.matching(variantStore.selectAll()).map { it.vcfId }.toSet()

        logger.info("Identifying PON filtered variants")
        val ponFiltered = ponStore.matching(variantStore.selectAll()).map { it.vcfId }.toSet()

        logger.info("Applying initial soft filters")
        val initialFilters = SoftFilterStore(config.filterConfig, variantStore.selectAll(), ponFiltered, hotspots)

        logger.info("Finding assembly links")
        val assemblyLinks: LinkStore = AssemblyLink(variantStore.selectAll())

        logger.info("Finding transitive links")
        val alternatePaths: Collection<AlternatePath> = AlternatePath(assemblyLinks, variantStore)
        val alternatePathsStringsByVcfId = alternatePaths.associate { x -> Pair(x.vcfId, x.pathString()) }
        val transitiveLinks = LinkStore(alternatePaths.flatMap { x -> x.transitiveLinks() })

        logger.info("Paired break end de-duplication")
        val dedupPair = DedupPair(initialFilters, alternatePaths, variantStore)
        val softFiltersAfterPairedDedup = initialFilters.update(dedupPair.duplicates, dedupPair.rescue)

        logger.info("Single break end de-duplication")
        val dedupSingle = DedupSingle(variantStore, softFiltersAfterPairedDedup)
        val softFiltersAfterSingleDedup = softFiltersAfterPairedDedup.update(dedupSingle.duplicates, setOf())

        logger.info("Finding double stranded break links")
        val dsbLinks = DsbLink(variantStore, assemblyLinks, softFiltersAfterSingleDedup.duplicates())
        val combinedLinks = LinkStore(assemblyLinks, transitiveLinks, dsbLinks)

        logger.info("Rescuing linked variants")
        val linkRescues = LinkRescue(combinedLinks, softFiltersAfterSingleDedup, variantStore)
        val finalFilters: SoftFilterStore = softFiltersAfterSingleDedup.update(setOf(), linkRescues.rescues)

        logger.info("Writing file: ${config.outputVcf}")
        fileWriter.writeHeader(fileReader.fileHeader)
        for (variant in variantStore.selectAll()) {

            val localLinkedBy = combinedLinks[variant.vcfId]
            val remoteLinkedBy = combinedLinks[variant.mateId]
            val altPath = alternatePathsStringsByVcfId[variant.vcfId]

            val filters = finalFilters.filters(variant.vcfId, variant.mateId)
            fileWriter.writeVariant(variant.context(localLinkedBy, remoteLinkedBy, altPath, hotspots.contains(variant.vcfId), filters))
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
                structuralVariants.add(structuralVariant.realign(refGenome))
            }
        }

        val mateIsValidOrNull = { x: StructuralVariantContext -> x.mateId?.let { unfiltered.contains(it) } != false }
        return structuralVariants.filter { x -> !hardFilter.contains(x.vcfId) && mateIsValidOrNull(x) }
    }


    override fun close() {
        refGenome.close()
        fileReader.close()
        fileWriter.close()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}

