package com.hartwig.hmftools.cider

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.cider.AsyncBamReader.processBam
import com.hartwig.hmftools.cider.VDJSequenceTsvWriter.writeVDJSequences
import com.hartwig.hmftools.cider.alignment.AlignmentAnnotation
import com.hartwig.hmftools.cider.alignment.AlignmentAnnotator
import com.hartwig.hmftools.cider.alignment.AlignmentStatus
import com.hartwig.hmftools.cider.genes.IgTcrConstantDiversityRegion
import com.hartwig.hmftools.cider.primer.PrimerTsvFile
import com.hartwig.hmftools.cider.primer.VdjPrimerMatch
import com.hartwig.hmftools.cider.primer.VdjPrimerMatchTsv
import com.hartwig.hmftools.cider.primer.VdjPrimerMatcher
import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.common.genome.region.GenomeRegions
import com.hartwig.hmftools.common.utils.config.ConfigBuilder
import com.hartwig.hmftools.common.utils.config.ConfigUtils
import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import com.hartwig.hmftools.common.utils.config.VersionInfo
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.SAMFileWriterFactory
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.cram.ref.ReferenceSource
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.IOException
import java.time.Duration
import java.time.Instant
import java.time.LocalDate
import java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME
import java.util.*
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import java.util.concurrent.Future
import kotlin.system.exitProcess

class CiderApplication(configBuilder: ConfigBuilder)
{
    private val mParams = CiderParams(configBuilder)

    @Throws(IOException::class, InterruptedException::class)
    fun run(args: Array<String>): Int
    {
        val runDate = LocalDate.now()

        val versionInfo = VersionInfo("cider.version")
        sLogger.info("Cider version: {}, build timestamp: {}",
            versionInfo.version(),
            versionInfo.buildTime().format(ISO_ZONED_DATE_TIME))

        sLogger.info("run date: {}", runDate)
        sLogger.info("run args: {}", args.joinToString(" "))

        FileWriterUtils.checkCreateOutputDir(mParams.outputDir)
        val start = Instant.now()

        val ciderGeneDatastore: ICiderGeneDatastore = CiderGeneDatastore(
            CiderGeneDataLoader.loadAnchorTemplates(mParams.refGenomeVersion),
            CiderGeneDataLoader.loadConstantDiversityRegions(mParams.refGenomeVersion))

        val candidateBlosumSearcher = AnchorBlosumSearcher(
            ciderGeneDatastore,
            CiderConstants.CANDIDATE_MIN_PARTIAL_ANCHOR_AA_LENGTH)

        val readProcessor = CiderReadScreener(
            ciderGeneDatastore,
            candidateBlosumSearcher,
            CiderConstants.MAX_READ_DISTANCE_FROM_ANCHOR,
            mParams.approxMaxFragmentLength)

        readBamFile(readProcessor, ciderGeneDatastore)
        writeCiderBam(readProcessor.allMatchedReads)

        val vjReadLayoutAdaptor = VJReadLayoutBuilder(mParams.numBasesToTrim, mParams.minBaseQuality)
        val layoutBuildResults = buildLayouts(vjReadLayoutAdaptor, readProcessor.vjReadCandidates, mParams.threadCount)

        val vdjBuilderBlosumSearcher = AnchorBlosumSearcher(
            ciderGeneDatastore,
            CiderConstants.VDJ_MIN_PARTIAL_ANCHOR_AA_LENGTH
        )

        val vdjSeqBuilder = VDJSequenceBuilder(
            vjReadLayoutAdaptor, vdjBuilderBlosumSearcher, mParams.minBaseQuality.toByte(),
            CiderConstants.MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES
        )

        val vdjSequences: List<VDJSequence> = vdjSeqBuilder.buildVDJSequences(layoutBuildResults.mapValues { (_, v) -> v.layouts })
        var primerMatchList: List<VdjPrimerMatch> = emptyList()

        if (mParams.primerCsv != null)
        {
            val primerList = PrimerTsvFile.load(mParams.primerCsv!!)
            // if we are provided a list of primers, match those against the input
            val vdjPrimerMatcher = VdjPrimerMatcher(mParams.primerMismatchMax)
            primerMatchList = vdjPrimerMatcher.matchVdjPrimer(vdjSequences, primerList)

            // write out the primer matches
            VdjPrimerMatchTsv.writePrimerMatches(mParams.outputDir, mParams.sampleId, primerMatchList)
        }

        val vdjAnnotator = VdjAnnotator(vjReadLayoutAdaptor, vdjBuilderBlosumSearcher)
        val alignmentAnnotations: Collection<AlignmentAnnotation>

        if (true)
        {
            // we need to filter out VDJ sequences that already match reference. In this version we avoid running alignment on those
            val filteredVdjs = vdjSequences.filter { vdj -> !vdjAnnotator.vdjMatchesRef(vdj) }

            // perform a GC collection before running alignment. This is to reduce memory used by JVM
            System.gc()

            val alignmentAnnotator = AlignmentAnnotator()
            alignmentAnnotations = alignmentAnnotator.runAnnotate(mParams.sampleId, filteredVdjs, mParams.outputDir, mParams.threadCount)
        }
        else
        {
            alignmentAnnotations = emptyList()
        }

        var vdjAnnotations: List<VdjAnnotation> = vdjAnnotator.sortAndAnnotateVdjs(vdjSequences, alignmentAnnotations, primerMatchList)

        // apply hard filters
        vdjAnnotations = vdjAnnotations.filter { vdjAnnotation -> passesHardFilters(vdjAnnotation) }

        writeVDJSequences(mParams.outputDir, mParams.sampleId, vdjAnnotations)

        // write the stats per locus
        CiderLocusStatsWriter.writeLocusStats(mParams.outputDir, mParams.sampleId, layoutBuildResults, vdjAnnotations)

        val finish: Instant = Instant.now()
        val seconds: Long = Duration.between(start, finish).seconds
        sLogger.info("CIDER run complete. Time taken: {}m {}s", seconds / 60, seconds % 60)
        return 0
    }

    @Throws(InterruptedException::class, IOException::class)
    fun readBamFile(readProcessor: CiderReadScreener, ciderGeneDatastore: ICiderGeneDatastore)
    {
        val readerFactory = readerFactory(mParams)
        val asyncBamRecordHander: (SAMRecord) -> Unit = { samRecord: SAMRecord ->
            readProcessor.asyncProcessSamRecord(samRecord)
        }

        val genomeRegions = TreeSet<GenomeRegion>()

        // first add all the VJ anchor locations
        for (anchorGenomeLoc: VJAnchorGenomeLocation in ciderGeneDatastore.getVjAnchorGeneLocations())
        {
            require(anchorGenomeLoc.genomeLocation.inPrimaryAssembly)
            genomeRegions.add(GenomeRegions.create(
                anchorGenomeLoc.chromosome,
                anchorGenomeLoc.start - mParams.approxMaxFragmentLength,
                anchorGenomeLoc.end + mParams.approxMaxFragmentLength))
        }

        // then add all the constant / diversity region genome locations
        for (region: IgTcrConstantDiversityRegion in ciderGeneDatastore.getIgConstantDiversityRegions())
        {
            require(region.genomeLocation.inPrimaryAssembly)
            genomeRegions.add(GenomeRegions.create(
                region.genomeLocation.chromosome,
                region.genomeLocation.posStart - mParams.approxMaxFragmentLength,
                region.genomeLocation.posEnd + mParams.approxMaxFragmentLength))
        }

        processBam(mParams.bamPath, readerFactory, genomeRegions, asyncBamRecordHander, mParams.threadCount)
        sLogger.info("found {} VJ read records", readProcessor.allMatchedReads.size)
    }

    // Build the consensus layout of all reads
    private fun buildLayouts(
        vjReadLayoutAdaptor: VJReadLayoutBuilder, readCandidates: Collection<VJReadCandidate>, threadCount: Int)
        : Map<VJGeneType, VJReadLayoutBuilder.LayoutBuildResult>
    {
        val geneTypes = VJGeneType.values()

        // use a EnumMap such that the keys are ordered by the declaration
        val layoutResults: MutableMap<VJGeneType, VJReadLayoutBuilder.LayoutBuildResult> = EnumMap(VJGeneType::class.java)

        val namedThreadFactory = ThreadFactoryBuilder().setNameFormat("worker-%d").build()
        val executorService = Executors.newFixedThreadPool(threadCount, namedThreadFactory)

        try
        {
            val futures: MutableMap<VJGeneType, Future<VJReadLayoutBuilder.LayoutBuildResult>> = EnumMap(VJGeneType::class.java)

            for (geneType in geneTypes)
            {
                val readsOfGeneType = readCandidates
                    .filter { o: VJReadCandidate -> o.vjGeneType === geneType }
                    .toList()

                val workerTask: Callable<VJReadLayoutBuilder.LayoutBuildResult> = Callable {
                    vjReadLayoutAdaptor.buildLayouts(
                        geneType, readsOfGeneType, CiderConstants.LAYOUT_MIN_READ_OVERLAP_BASES,
                        mParams.maxReadCountPerGene, mParams.maxLowQualBaseFraction)
                }
                futures[geneType] = executorService.submit(workerTask)
            }

            for ((geneType, future) in futures)
            {
                layoutResults[geneType] = future.get()
            }
        }
        finally
        {
            // we must do this to make sure application will exit on exception
            executorService.shutdown()
        }

        // give each an ID
        var nextId = 1
        for ((_, layoutResult) in layoutResults)
        {
            for (layout in layoutResult.layouts)
            {
                layout.id = (nextId++).toString()
            }
        }

        // write all the layouts
        VJReadLayoutFile.writeLayouts(mParams.outputDir, mParams.sampleId, layoutResults.mapValues { (_, v) -> v.layouts })
        return layoutResults
    }

    @Throws(IOException::class)
    fun writeCiderBam(samRecords: Collection<SAMRecord?>)
    {
        if (!mParams.writeFilteredBam) return
        val outBamPath = mParams.outputDir + "/" + mParams.sampleId + ".cider.bam"
        val readerFactory = readerFactory(mParams)
        var samFileHeader: SAMFileHeader
        readerFactory.open(File(mParams.bamPath)).use { samReader -> samFileHeader = samReader.fileHeader }
        SAMFileWriterFactory().makeBAMWriter(
            samFileHeader, false, File(outBamPath)
        ).use { bamFileWriter ->
            for (r in samRecords)
            {
                bamFileWriter.addAlignment(r)
            }
        }
    }

    fun passesHardFilters(vdjAnnotation: VdjAnnotation) : Boolean
    {
        if (!mParams.reportMatchRefSeq && vdjAnnotation.filters.contains(VdjAnnotation.Filter.MATCHES_REF))
        {
            return false
        }

        // filter out sequences that are too short and only has V or J
        if (vdjAnnotation.filters.contains(VdjAnnotation.Filter.MIN_LENGTH))
        {
            if (vdjAnnotation.alignmentAnnotation == null)
            {
                // if alignment is not run, we use NO_V_ANCHOR / NO_J_ANCHOR
                if (vdjAnnotation.filters.contains(VdjAnnotation.Filter.NO_V_ANCHOR) ||
                    vdjAnnotation.filters.contains(VdjAnnotation.Filter.NO_J_ANCHOR))
                {
                    return false
                }
            }
            else if (vdjAnnotation.alignmentAnnotation!!.alignmentStatus == AlignmentStatus.V_ONLY ||
                    vdjAnnotation.alignmentAnnotation!!.alignmentStatus == AlignmentStatus.J_ONLY)
            {
                return false
            }
        }

        return true
    }

    companion object
    {
        val sLogger = LogManager.getLogger(CiderApplication::class.java)
        private fun readerFactory(params: CiderParams): SamReaderFactory
        {
            val readerFactory = SamReaderFactory.make()
            return if (params.refGenomePath != null)
            {
                readerFactory.referenceSource(ReferenceSource(File(params.refGenomePath!!)))
            } else readerFactory
        }

        @Throws(IOException::class, InterruptedException::class)
        @JvmStatic
        fun main(args: Array<String>)
        {
            val configBuilder = ConfigBuilder("Cider")
            CiderParams.registerConfig(configBuilder)
            ConfigUtils.addLoggingOptions(configBuilder)
            configBuilder.checkAndParseCommandLine(args)
            val ciderApplication = CiderApplication(configBuilder)
            exitProcess(ciderApplication.run(args))
        }
    }
}