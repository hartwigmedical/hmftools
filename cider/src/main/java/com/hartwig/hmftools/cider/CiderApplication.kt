package com.hartwig.hmftools.cider

import com.beust.jcommander.JCommander
import com.beust.jcommander.ParameterException
import com.beust.jcommander.ParametersDelegate
import com.beust.jcommander.UnixStyleUsageFormatter
import com.hartwig.hmftools.cider.AsyncBamReader.processBam
import com.hartwig.hmftools.cider.VDJSequenceTsvWriter.writeVDJSequences
import com.hartwig.hmftools.cider.VJReadLayoutFile.writeLayouts
import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.cider.primer.*
import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.common.genome.region.GenomeRegions
import com.hartwig.hmftools.common.utils.FileWriterUtils
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import com.hartwig.hmftools.common.utils.config.LoggingOptions
import com.hartwig.hmftools.common.utils.version.VersionInfo
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
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter.ISO_ZONED_DATE_TIME
import java.util.*
import kotlin.collections.ArrayList

class CiderApplication
{
    // add the options
    @ParametersDelegate
    private val mParams = CiderParams()

    // add to the logging options
    @ParametersDelegate
    private val mLoggingOptions = LoggingOptions()

    @Throws(IOException::class, InterruptedException::class)
    fun run(): Int
    {
        mLoggingOptions.setLogLevel()
        val versionInfo = VersionInfo("cider.version")
        sLogger.info("Cider version: {}, build timestamp: {}",
            versionInfo.version(),
            versionInfo.buildTime()!!.format(ISO_ZONED_DATE_TIME))

        if (!mParams.isValid)
        {
            sLogger.error(" invalid config, exiting")
            return 1
        }
        FileWriterUtils.checkCreateOutputDir(mParams.outputDir)
        val start = Instant.now()

        val ciderGeneDatastore: ICiderGeneDatastore = CiderGeneDatastore(
            CiderGeneDataLoader.loadAnchorTemplateTsv(mParams.refGenomeVersion),
            CiderGeneDataLoader.loadConstantRegionGenes(mParams.refGenomeVersion, mParams.ensemblDataDir))

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
        val layoutMap = buildLayouts(vjReadLayoutAdaptor, readProcessor.vjReadCandidates)

        val vdjBuilderBlosumSearcher = AnchorBlosumSearcher(
            ciderGeneDatastore,
            CiderConstants.VDJ_MIN_PARTIAL_ANCHOR_AA_LENGTH
        )

        val vdjSeqBuilder = VDJSequenceBuilder(
            vjReadLayoutAdaptor, vdjBuilderBlosumSearcher, mParams.minBaseQuality.toByte(),
            CiderConstants.MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES
        )

        val vdjSequences: List<VDJSequence> = vdjSeqBuilder.buildVDJSequences(layoutMap)
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
        val vdjAnnotations: List<VdjAnnotation> = vdjAnnotator.sortAndAnnotateVdjs(vdjSequences, primerMatchList)

        writeVDJSequences(mParams.outputDir, mParams.sampleId, vdjAnnotations, mParams.reportMatchRefSeq, true)

        val finish = Instant.now()
        val seconds = Duration.between(start, finish).seconds
        sLogger.info("CIDER run complete, time taken: {}m {}s", seconds / 60, seconds % 60)
        return 0
    }

    @Throws(InterruptedException::class, IOException::class)
    fun readBamFile(readProcessor: CiderReadScreener, ciderGeneDatastore: ICiderGeneDatastore)
    {
        val readerFactory = readerFactory(mParams)
        val asyncBamRecordHander: (SAMRecord) -> Unit = { samRecord: SAMRecord ->
            readProcessor.asyncProcessSamRecord(samRecord)
        }

        val genomeRegions = ArrayList<GenomeRegion>()

        // first add all the VJ anchor locations
        for (anchorGenomeLoc: VJAnchorGenomeLocation in ciderGeneDatastore.getVjAnchorGeneLocations())
        {
            genomeRegions.add(GenomeRegions.create(
                anchorGenomeLoc.chromosome,
                anchorGenomeLoc.start - mParams.approxMaxFragmentLength,
                anchorGenomeLoc.end + mParams.approxMaxFragmentLength))
        }

        // then add all the constant region genome locations
        for (constantRegion: IgTcrConstantRegion in ciderGeneDatastore.getIgConstantRegions())
        {
            genomeRegions.add(GenomeRegions.create(
                constantRegion.genomeLocation.chromosome,
                constantRegion.genomeLocation.posStart - mParams.approxMaxFragmentLength,
                constantRegion.genomeLocation.posEnd + mParams.approxMaxFragmentLength))
        }
        processBam(mParams.bamPath, readerFactory, genomeRegions, asyncBamRecordHander, mParams.threadCount)
        sLogger.info("found {} VJ read records", readProcessor.allMatchedReads.size)
    }

    private fun buildLayouts(
        vjReadLayoutAdaptor: VJReadLayoutBuilder,
        readCandidates: Collection<VJReadCandidate>
    ): Map<VJGeneType, List<ReadLayout>>
    {
        // now build the consensus overlay sequences
        //var geneTypes = new VJGeneType[] { VJGeneType.IGHV, VJGeneType.IGHJ };
        val geneTypes = VJGeneType.values()

        // use a EnumMap such that the keys are ordered by the declaration
        val layoutMap: MutableMap<VJGeneType, List<ReadLayout>> = EnumMap(VJGeneType::class.java)

        for (geneType in geneTypes)
        {
            val readsOfGeneType = readCandidates
                .filter({ o: VJReadCandidate -> o.vjGeneType === geneType })
                .toList()

            val readLayouts = vjReadLayoutAdaptor.buildLayouts(
                geneType, readsOfGeneType, CiderConstants.LAYOUT_MIN_READ_OVERLAP_BASES)

            layoutMap[geneType] = readLayouts.sortedByDescending({ layout: ReadLayout -> layout.reads.size })
        }

        // give each an ID
        var nextId = 1
        for ((_, layoutList) in layoutMap)
        {
            for (layout in layoutList)
            {
                layout.id = (nextId++).toString()
            }
        }
        writeLayouts(mParams.outputDir, mParams.sampleId, layoutMap)
        return layoutMap
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
            sLogger.info("{}", LocalDateTime.now())
            sLogger.info("args: {}", java.lang.String.join(" ", *args))
            val ciderApplication = CiderApplication()
            val commander = JCommander.newBuilder()
                .addObject(ciderApplication)
                .build()

            // use unix style formatter
            commander.usageFormatter = UnixStyleUsageFormatter(commander)
            // help message show in order parameters are declared
            commander.parameterDescriptionComparator = DeclaredOrderParameterComparator(CiderApplication::class.java)
            try
            {
                commander.parse(*args)
            } catch (e: ParameterException)
            {
                println("Unable to parse args: " + e.message)
                commander.usage()
                System.exit(1)
            }

            // set all thread exception handler
            Thread.setDefaultUncaughtExceptionHandler(
                { t: Thread, e: Throwable ->
                    sLogger.error("[{}]: uncaught exception: {}", t, e)
                    e.printStackTrace(System.err)
                    System.exit(1)
                })

            System.exit(ciderApplication.run())
        }
    }
}