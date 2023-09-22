package com.hartwig.hmftools.cider.blastn

import com.beust.jcommander.JCommander
import com.beust.jcommander.ParameterException
import com.beust.jcommander.ParametersDelegate
import com.beust.jcommander.UnixStyleUsageFormatter
import com.google.common.collect.Multimap
import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.cider.VDJSequenceTsvWriter
import com.hartwig.hmftools.cider.VdjAnnotation
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import com.hartwig.hmftools.common.utils.config.LoggingOptions
import com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader
import com.hartwig.hmftools.common.utils.version.VersionInfo
import htsjdk.samtools.util.SequenceUtil
import org.apache.commons.csv.CSVFormat
import org.apache.logging.log4j.LogManager
import java.io.IOException
import java.time.LocalDate
import java.time.format.DateTimeFormatter
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import java.util.concurrent.Future

// find all of the sequences in the reference genome that matches with
// the cider sequences. And write them all to a file
class CiderRefGenomeExtractor
{
    // add the options
    @ParametersDelegate
    private val mParams = RefGenomeExtractorParams()

    // add to the logging options
    @ParametersDelegate
    private val mLoggingOptions = LoggingOptions()

    @Throws(IOException::class, InterruptedException::class)
    fun run(args: Array<String>): Int
    {
        val runDate = LocalDate.now()

        mLoggingOptions.setLogLevel()
        val versionInfo = VersionInfo("cider.version")
        sLogger.info(
            "Cider version: {}, build timestamp: {}",
            versionInfo.version(),
            versionInfo.buildTime().format(DateTimeFormatter.ISO_ZONED_DATE_TIME)
        )

        val vdjTsvs = ArrayList<String>()

        createBufferedReader(mParams.config).use { reader ->
            var line: String?

            while (true)
            {
                line = reader.readLine() ?: break
                vdjTsvs.add(line)
            }
        }

        val refGenomeRegionCollator = RefGenomeRegionCollator()

        if (mParams.inputRefGenomeRegions != null)
        {
            refGenomeRegionCollator.readFromTsv(mParams.inputRefGenomeRegions!!)
        }

        // we want to multi thread the results
        val namedThreadFactory = ThreadFactoryBuilder().setNameFormat("worker-%d").build()
        val executorService = Executors.newFixedThreadPool(mParams.threadCount, namedThreadFactory)

        try
        {
            val futures = ArrayList<Future<Unit>>()

            // for each vdj tsv, we run blastn on it
            for (vdjTsv in vdjTsvs)
            {
                val task = Callable { processVdjTsv(vdjTsv, refGenomeRegionCollator) }
                futures.add(executorService.submit(task))
            }

            var i = 0
            for (fut in futures)
            {
                fut.get()

                if (++i == 10)
                {
                    synchronized(refGenomeRegionCollator) {
                        refGenomeRegionCollator.writeToTsv(mParams.outputRefGenomeRegions)
                    }
                    i = 0
                }
            }
        }
        finally
        {
            // we must do this to make sure application will exit on exception
            executorService.shutdown()
        }

        // write out the ref genome regions found
        refGenomeRegionCollator.writeToTsv(mParams.outputRefGenomeRegions)

        return 0
    }

    fun processVdjTsv(vdjTsv: String, refGenomeRegionCollator: RefGenomeRegionCollator)
    {
        var key = 0

        val vdjSequences = HashMap<Int, String>()

        createBufferedReader(vdjTsv).use { reader ->
            val format = CSVFormat.Builder.create()
                .setDelimiter('\t')
                .setRecordSeparator('\n')
                .setHeader().setSkipHeaderRecord(true)
                .build()

            for (record in format.parse(reader))
            {
                val fullSeq = record[VDJSequenceTsvWriter.Column.fullSeq.toString()]
                val filter = record[VDJSequenceTsvWriter.Column.filter.toString()]

                if (filter.contains(VdjAnnotation.Filter.MATCHES_REF.toString()) ||
                    filter.contains(VdjAnnotation.Filter.DUPLICATE.toString()) )
                {
                    continue
                }

                // we use full sequence to do blastn
                vdjSequences[key++] = fullSeq
            }
        }

        // parse out the sample id
        var sampleId = vdjTsv
        var i: Int = Math.max(sampleId.lastIndexOf('/'), sampleId.lastIndexOf('\\'))
        if (i != -1)
        {
            sampleId = sampleId.substring(i + 1)
        }
        i = sampleId.indexOf('.')
        if (i != -1)
        {
            sampleId = sampleId.substring(0, i)
        }

        // run blastn
        val blastnMatches: Multimap<Int, BlastnMatch> = BlastnRunner.runBlastn(
            sampleId, mParams.blast, mParams.blastDb, vdjSequences, mParams.tempDir, 1,
            expectedValueCutoff = 0.1)

        synchronized(refGenomeRegionCollator) {

            // add all blastn matches
            for ((key, blastnMatch) in blastnMatches.entries())
            {
                // add each match into the list
                // data class RefGenomeRegion(val contig: String, var start: Int, var endExclusive: Int, var sequence: String)
                val refGenomeRegion: RefGenomeRegion = if (blastnMatch.subjectFrame == Strand.FORWARD)
                {
                    RefGenomeRegion(
                        contig = blastnMatch.subjectTitle,
                        start = blastnMatch.subjectAlignStart,
                        endExclusive = blastnMatch.subjectAlignEnd + 1,
                        sequence = blastnMatch.alignedPartOfSubjectSeq.replace("-", "")
                    )
                } else
                {
                    RefGenomeRegion(
                        contig = blastnMatch.subjectTitle,
                        start = blastnMatch.subjectAlignEnd,
                        endExclusive = blastnMatch.subjectAlignStart + 1,
                        sequence = SequenceUtil.reverseComplement(blastnMatch.alignedPartOfSubjectSeq.replace("-", ""))
                    )
                }

                refGenomeRegionCollator.addRefGenomeRegion(refGenomeRegion)

                //sLogger.debug("added ref genome region: {}", refGenomeRegion)
            }
        }
    }

    companion object
    {
        val sLogger = LogManager.getLogger(CiderRefGenomeExtractor::class.java)

        @Throws(IOException::class, InterruptedException::class)
        @JvmStatic
        fun main(args: Array<String>)
        {
            val refGenomeExtractor = CiderRefGenomeExtractor()
            val commander = JCommander.newBuilder()
                .addObject(refGenomeExtractor)
                .build()

            // use unix style formatter
            commander.usageFormatter = UnixStyleUsageFormatter(commander)
            // help message show in order parameters are declared
            commander.parameterDescriptionComparator = DeclaredOrderParameterComparator(CiderRefGenomeExtractor::class.java)
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

            System.exit(refGenomeExtractor.run(args))
        }
    }
}