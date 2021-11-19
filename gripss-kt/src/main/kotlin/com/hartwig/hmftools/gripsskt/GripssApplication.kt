package com.hartwig.hmftools.gripsskt

import com.hartwig.hmftools.common.utils.version.VersionInfo
import com.hartwig.hmftools.gripsskt.GripssApplication.Companion.logger
import com.hartwig.hmftools.gripsskt.dedup.DedupPair
import com.hartwig.hmftools.gripsskt.dedup.DedupSingle
import com.hartwig.hmftools.gripsskt.link.AlternatePath
import com.hartwig.hmftools.gripsskt.link.AssemblyLink
import com.hartwig.hmftools.gripsskt.link.DsbLink
import com.hartwig.hmftools.gripsskt.link.LinkRescue
import com.hartwig.hmftools.gripsskt.store.*
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.variant.variantcontext.VariantContextComparator
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeader
import org.apache.commons.cli.*
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.IOException
import kotlin.system.exitProcess

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
        exitProcess(1)
    } catch (e: ParseException) {
        logger.warn(e)
        val formatter = HelpFormatter()
        formatter.printHelp("gripss", GripssConfig.createHelpOptions())
        exitProcess(1)
    }
}

class GripssApplication(private val config: GripssConfig) : AutoCloseable, Runnable {
    companion object {
        val logger = LogManager.getLogger(this::class.java)
        val version = VersionInfo("gripss.version")
    }

    private val startTime = System.currentTimeMillis()
    private val fileReader = VCFFileReader(File(config.inputVcf), false)
    private val dictionary = fileReader.fileHeader.sequenceDictionary
    private val fileWriter = GripssVCF(config.outputVcf, dictionary)
    private val refGenome = IndexedFastaSequenceFile(File(config.refGenome))

    override fun run() {
        val inputHeader = fileReader.fileHeader!!
        assertInfoLine(inputHeader, BEALN)
        assertInfoLine(inputHeader, REPEAT_MASKER_REPEAT_CLASS)
        assertInfoLine(inputHeader, REPEAT_MASKER_REPEAT_TYPE)

        logger.info("Config ${config.filterConfig}")
        val contigComparator = ContigComparator(dictionary)

        val inputSampleNames = inputHeader.genotypeSamples!!
        val inputSampleOrdinals = config.sampleOrdinals(inputSampleNames)
        val outputSampleNames = mutableListOf<String>()

        if (inputSampleOrdinals.first != -1) {
            val referenceSample = inputSampleNames[inputSampleOrdinals.first]
            logger.info("Using $referenceSample as reference sample")
            outputSampleNames.add(referenceSample)
        } else {
            logger.info("Running in tumor-only mode")
        }
        val tumorSample = inputSampleNames[inputSampleOrdinals.second]
        logger.info("Using $tumorSample as tumor sample")
        outputSampleNames.add(tumorSample)

        logger.info("Reading hotspot file: ${config.pairedHotspotFile}")
        val pairedHotspots = Breakpoint.fromBedpeFile(config.pairedHotspotFile, contigComparator)
        val hotspotStore = HotspotStore(contigComparator, pairedHotspots)
        val hotspotRescue = hotspotRescue(hotspotStore)

        logger.info("Reading VCF file: ${config.inputVcf}")
        val variantStore = VariantStore(hardFilterAndRealign(fileReader, inputSampleOrdinals, hotspotRescue, contigComparator))
        val hotspots = variantStore.selectAll().filter(hotspotRescue).map { x -> x.vcfId }.toSet()

        logger.info("Reading PON files: ${config.singlePonFile} ${config.pairedPonFile}")
        val ponFiltered = ponFiltered(contigComparator, variantStore.selectAll())

        logger.info("Applying initial soft filters")
        val initialFilters = SoftFilterStore(config.filterConfig, contigComparator, variantStore.selectAll(), ponFiltered, hotspots)

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
        fileWriter.writeHeader(version.version(), fileReader.fileHeader, outputSampleNames)
        for (variant in variantStore.selectAll()) {

            val localLinkedBy = combinedLinks[variant.vcfId]
            val remoteLinkedBy = combinedLinks[variant.mateId]
            val altPath = alternatePathsStringsByVcfId[variant.vcfId]

            val filters = finalFilters.filters(variant.vcfId, variant.mateId)
            fileWriter.writeVariant(variant.context(localLinkedBy, remoteLinkedBy, altPath, hotspots.contains(variant.vcfId), filters))
        }
    }

    private fun hotspotRescue(hotspotStore: HotspotStore): (StructuralVariantContext) -> Boolean {
        val appropriateSoftFilters = { x: StructuralVariantContext -> !x.polyATHomologyFilter() }
        val minDistanceFilter = { x: StructuralVariantContext -> !x.isTooShortToRescue }
        return { variant -> minDistanceFilter(variant) && appropriateSoftFilters(variant) && hotspotStore.contains(variant) }
    }

    private fun ponFiltered(contigComparator: ContigComparator, variants: List<StructuralVariantContext>): Set<String> {
        val breakends = Breakend.fromBedFile(config.singlePonFile)
        val breakpoints = Breakpoint.fromBedpeFile(config.pairedPonFile, contigComparator)
        val ponStore = LocationStore(contigComparator, breakends, breakpoints, config.filterConfig.ponDistance)

        logger.info("Applying PON file")
        return variants.filter { ponStore.contains(it) }.map { it.vcfId }.toSet()
    }

    private fun hardFilterAndRealign(fileReader: VCFFileReader, ordinals: Pair<Int, Int>, hotspotFilter: (StructuralVariantContext) -> Boolean, contigComparator: ContigComparator): List<StructuralVariantContext> {
        val hardFilter: MutableSet<String> = mutableSetOf()
        val validVariantsById: MutableMap<String, StructuralVariantContext> = mutableMapOf()

        // Read in all non hard filtered
        for (variantContext in fileReader) {
            val structuralVariant = StructuralVariantContext(variantContext, ordinals.first, ordinals.second)
            val isHotspot = hotspotFilter(structuralVariant)
            val isMateHardFiltered = hardFilter.contains(structuralVariant.vcfId)
            val isHardFiltered = isMateHardFiltered || structuralVariant.isHardFilter(config.filterConfig, contigComparator, isHotspot)

            if (isHardFiltered) {
                structuralVariant.mateId?.let {
                    hardFilter.add(it)
                    validVariantsById.remove(it)
                }
            } else {
                validVariantsById[variantContext.id] = structuralVariant
            }
        }

        // Realignment
        val vcfIds = validVariantsById.keys
        for (vcfId in vcfIds) {
            val original = validVariantsById[vcfId]!!
            val realigned = original.realign(refGenome, contigComparator)
            if (realigned.context.realigned()) {
                validVariantsById[vcfId] = realigned
                realigned.mateId?.let { validVariantsById[it] }?.let { x ->
                    validVariantsById[x.vcfId] = x.realignRemote(realigned)
                }
            }
        }

        val contextComparator = VariantContextComparator(dictionary)
        val comparator: Comparator<StructuralVariantContext> = Comparator { x, y -> contextComparator.compare(x.context, y.context) }
        val mateIsValidOrNull = { x: StructuralVariantContext -> x.mateId?.let { validVariantsById.contains(it) } != false }
        return validVariantsById.values.filter { x -> !hardFilter.contains(x.vcfId) && mateIsValidOrNull(x) }.sortedWith(comparator)
    }

    private fun assertInfoLine(header: VCFHeader, tag: String) {
        if (!header.hasInfoLine(tag)) {
            throw IllegalArgumentException("Supplied VCF must be enriched with $tag annotation")
        }
    }

    override fun close() {
        refGenome.close()
        fileReader.close()
        fileWriter.close()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }
}


