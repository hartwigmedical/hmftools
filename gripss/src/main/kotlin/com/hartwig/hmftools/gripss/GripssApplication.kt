package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.bedpe.Breakend
import com.hartwig.hmftools.bedpe.Breakpoint
import com.hartwig.hmftools.common.utils.version.VersionInfo
import com.hartwig.hmftools.gripss.GripssApplication.Companion.logger
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
import htsjdk.variant.variantcontext.VariantContextComparator
import htsjdk.variant.vcf.VCFFileReader
import org.apache.commons.cli.*
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.IOException

fun main(args: Array<String>) {

    @Throws(ParseException::class)
    fun createCommandLine(args: Array<String>, options: Options): CommandLine {
        val parser: CommandLineParser = DefaultParser()
        return parser.parse(options, args)
    }

    val options = GripssConfig.createOptions()
    try {
        val cmd = createCommandLine(args, options)
        val config = GripssConfig.createConfig(cmd)
        GripssApplication(config).use { x -> x.run() }
    } catch (e: IOException) {
        logger.warn(e)
    } catch (e: ParseException) {
        logger.warn(e)
        val formatter = HelpFormatter()
        formatter.printHelp("gripss", options)
    }
}

class GripssApplication(private val config: GripssConfig) : AutoCloseable, Runnable {
    companion object {
        const val PON_ADDITIONAL_DISTANCE = 1

        val logger = LogManager.getLogger(this::class.java)
        val version = VersionInfo("gripss.version")
    }

    private val startTime = System.currentTimeMillis()
    private val fileReader = VCFFileReader(File(config.inputVcf), false)
    private val dictionary = fileReader.fileHeader.sequenceDictionary
    private val fileWriter = GripssVCF(config.outputVcf, dictionary)
    private val refGenome = IndexedFastaSequenceFile(File(config.refGenome))

    override fun run() {
        logger.info("Config ${config.filterConfig}")
        val contigComparator = ContigComparator(dictionary)

        val sampleNames =  fileReader.fileHeader!!.genotypeSamples!!
        val sampleOrdinals = sampleOrdinals(sampleNames);
        logger.info("Using ${sampleNames[sampleOrdinals.first]} as reference, ${sampleNames[sampleOrdinals.second]} as tumor")

        logger.info("Reading hotspot file: ${config.pairedHotspotFile}")
        val hotspotStore = LocationStore(contigComparator, listOf(), Breakpoint.fromBedpeFile(config.pairedHotspotFile, contigComparator))
        val hotspotFilter = hotspotFilter(hotspotStore)

        logger.info("Reading VCF file: ${config.inputVcf}")
        val variantStore = VariantStore(hardFilterAndRealign(fileReader, sampleOrdinals, hotspotFilter, contigComparator))
        val hotspots = variantStore.selectAll().filter(hotspotFilter).map { x -> x.vcfId }.toSet()

        logger.info("Reading PON files: ${config.singlePonFile} ${config.pairedPonFile}")
        val ponFiltered = ponFiltered(contigComparator, variantStore.selectAll())

        logger.info("Applying initial soft filters")
        val initialFilters = SoftFilterStore(config.filterConfig, variantStore.selectAll(), ponFiltered, hotspots)

        logger.info("Finding assembly links")
        val assemblyLinks: LinkStore = AssemblyLink(variantStore.selectAll())

        logger.info("Finding transitive links")
        val alternatePaths: Collection<AlternatePath> = AlternatePath(assemblyLinks, variantStore)
        val alternatePathsStringsByVcfId = alternatePaths.associate { x -> Pair(x.vcfId, x.pathString()) }
        val transitiveLinks = LinkStore(alternatePaths.flatMap { x -> x.transitiveLinks() })
        val combinedTransitiveAssemblyLinks = LinkStore(assemblyLinks, transitiveLinks)

        logger.info("Paired break end de-duplication")
        val dedupPair = DedupPair(initialFilters, alternatePaths, variantStore)
        val softFiltersAfterPairedDedup = initialFilters.update(dedupPair.duplicates, dedupPair.rescue)

        logger.info("Single break end de-duplication")
        val dedupSingle = DedupSingle(variantStore, softFiltersAfterPairedDedup, combinedTransitiveAssemblyLinks)
        val softFiltersAfterSingleDedup = softFiltersAfterPairedDedup.update(dedupSingle.duplicates, setOf())

        logger.info("Finding double stranded break links")
        val dsbLinks = DsbLink(variantStore, assemblyLinks, softFiltersAfterSingleDedup.duplicates())

        logger.info("Rescuing linked variants")
        val dsbRescues = LinkRescue.rescueDsb(dsbLinks, softFiltersAfterSingleDedup, variantStore).rescues
        val dsbRescueMobileElements = LinkRescue.rescueDsbMobileElementInsertion(config.filterConfig, dsbLinks, softFiltersAfterSingleDedup, variantStore).rescues
        val assemblyRescues = LinkRescue.rescueAssembly(assemblyLinks, softFiltersAfterSingleDedup, variantStore).rescues
        val transitiveRescues = LinkRescue.rescueTransitive(transitiveLinks, softFiltersAfterSingleDedup, variantStore).rescues
        val allRescues = dsbRescues + dsbRescueMobileElements + assemblyRescues + transitiveRescues


        logger.info("Writing file: ${config.outputVcf}")
        val combinedLinks = LinkStore(combinedTransitiveAssemblyLinks, dsbLinks)
        val finalFilters: SoftFilterStore = softFiltersAfterSingleDedup.update(setOf(), allRescues)
        fileWriter.writeHeader(version.version(), fileReader.fileHeader)
        for (variant in variantStore.selectAll()) {

            val localLinkedBy = combinedLinks[variant.vcfId]
            val remoteLinkedBy = combinedLinks[variant.mateId]
            val altPath = alternatePathsStringsByVcfId[variant.vcfId]

            val filters = finalFilters.filters(variant.vcfId, variant.mateId)
            fileWriter.writeVariant(variant.context(localLinkedBy, remoteLinkedBy, altPath, hotspots.contains(variant.vcfId), filters))
        }
    }

    private fun sampleOrdinals(sampleNames: List<String>): Pair<Int, Int> {
        val normalOrdinal = if (config.reference.isEmpty()) 0 else sampleNames.indexOf(config.reference)
        if (normalOrdinal < 0) {
            throw IllegalArgumentException("Unable to locate sample ${config.reference} in supplied VCF")
        }
        val tumorOrdinal = if (config.tumor.isEmpty()) 1 else sampleNames.indexOf(config.tumor)
        if (tumorOrdinal < 0) {
            throw IllegalArgumentException("Unable to locate sample ${config.tumor} in supplied VCF")
        }

        if (tumorOrdinal == normalOrdinal) {
            throw ParseException("Tumor and reference must be different")
        }
        return Pair(normalOrdinal, tumorOrdinal)
    }

    private fun hotspotFilter(hotspotStore: LocationStore): (StructuralVariantContext) -> Boolean {
        val appropriateSoftFilters = { x: StructuralVariantContext -> !x.polyATHomologyFilter() }
        val minDistanceFilter = { x: StructuralVariantContext -> !x.isShort }
        return { variant -> minDistanceFilter(variant) && appropriateSoftFilters(variant) && hotspotStore.contains(variant) }
    }

    private fun ponFiltered(contigComparator: ContigComparator, variants: List<StructuralVariantContext>): Set<String> {
        val breakends = Breakend.fromBedFile(config.singlePonFile)
        val breakpoints = Breakpoint.fromBedpeFile(config.pairedPonFile, contigComparator)
        val ponStore = LocationStore(contigComparator, breakends, breakpoints, PON_ADDITIONAL_DISTANCE)

        logger.info("Applying PON file")
        return variants.filter { ponStore.contains(it) }.map { it.vcfId }.toSet()
    }

    private fun hardFilterAndRealign(fileReader: VCFFileReader, ordinals: Pair<Int, Int>, hotspotFilter: (StructuralVariantContext) -> Boolean, contigComparator: ContigComparator): List<StructuralVariantContext> {
        val unfiltered: MutableSet<String> = mutableSetOf()
        val hardFilter: MutableSet<String> = mutableSetOf()
        val structuralVariants: MutableList<StructuralVariantContext> = mutableListOf()

        for (variantContext in fileReader) {
            val structuralVariant = StructuralVariantContext(variantContext, ordinals.first, ordinals.second)
            val isHotspot = hotspotFilter(structuralVariant)
            val isMateFiltered = hardFilter.contains(structuralVariant.vcfId)
            val isHardFiltered = isMateFiltered || (!isHotspot && structuralVariant.isHardFilter(config.filterConfig))

            if (isHardFiltered) {
                structuralVariant.mateId?.let { hardFilter.add(it) }
            } else {
                unfiltered.add(variantContext.id)
                structuralVariants.add(structuralVariant.realign(refGenome, contigComparator))
            }
        }

        val contextComparator = VariantContextComparator(dictionary)
        val comparator: Comparator<StructuralVariantContext> = Comparator { x, y -> contextComparator.compare(x.context, y.context) }
        val mateIsValidOrNull = { x: StructuralVariantContext -> x.mateId?.let { unfiltered.contains(it) } != false }
        return structuralVariants.filter { x -> !hardFilter.contains(x.vcfId) && mateIsValidOrNull(x) }.sortedWith(comparator)
    }


    override fun close() {
        refGenome.close()
        fileReader.close()
        fileWriter.close()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}


