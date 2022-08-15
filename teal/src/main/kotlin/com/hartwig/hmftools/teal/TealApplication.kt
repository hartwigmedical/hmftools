package com.hartwig.hmftools.teal

import com.beust.jcommander.*
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCountFile
import com.hartwig.hmftools.common.genome.gc.ImmutableGCBucket
import com.hartwig.hmftools.common.metrics.WGSMetricsFile
import com.hartwig.hmftools.common.purple.PurityContextFile
import com.hartwig.hmftools.common.utils.FileWriterUtils
import com.hartwig.hmftools.common.utils.config.LoggingOptions
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import com.hartwig.hmftools.common.utils.version.VersionInfo
import com.hartwig.hmftools.teal.breakend.BreakEndApp
import com.hartwig.hmftools.teal.telbam.TelbamApp
import com.hartwig.hmftools.teal.tellength.SampleType
import com.hartwig.hmftools.teal.tellength.TelLengthApp
import org.apache.logging.log4j.LogManager
import java.time.Duration
import java.time.Instant
import kotlin.system.exitProcess

class TealApplication
{
    class StandaloneMode
    {
        @ParametersDelegate
        var params = TealParams()

        @ParametersDelegate
        val loggingOptions = LoggingOptions()

        fun run(): Int
        {
            params.validate()
            loggingOptions.setLogLevel()

            if (!FileWriterUtils.checkCreateOutputDir(params.commonParams.outputDir))
            {
                logger.error("failed to create output directory({})", params.commonParams.outputDir)
                return 1
            }

            val versionInfo = VersionInfo("teal.version")
            logger.info("Teal version: {}", versionInfo.version())
            logger.info("starting telomeric analysis")
            logger.info("{}", params)
            val start = Instant.now()

            val tumorOnly = params.commonParams.tumorOnly()
            val germlineOnly = params.commonParams.referenceOnly()

            try
            {
                logger.info("creating telbam files")

                // first we generate the telbam files
                // we pass the params to the telbam app

                if (!tumorOnly)
                {
                    val germlineTelbamApp = TelbamApp()
                    germlineTelbamApp.params.bamFile = params.commonParams.referenceBamFile!!
                    germlineTelbamApp.params.refGenomeFile = params.commonParams.refGenomeFile
                    germlineTelbamApp.params.telbamFile = germlineTelbamPath()
                    germlineTelbamApp.params.tsvFile = germlineTelReadTsvPath()
                    germlineTelbamApp.params.threadCount = params.commonParams.threadCount
                    germlineTelbamApp.processBam()
                }

                if (!germlineOnly)
                {
                    val tumorTelbamApp = TelbamApp()
                    tumorTelbamApp.params.bamFile = params.commonParams.tumorBamFile!!
                    tumorTelbamApp.params.refGenomeFile = params.commonParams.refGenomeFile
                    tumorTelbamApp.params.telbamFile = tumorTelbamPath()
                    tumorTelbamApp.params.tsvFile = tumorTelReadTsvPath()
                    tumorTelbamApp.params.threadCount = params.commonParams.threadCount
                    tumorTelbamApp.processBam()
                }

                var germlineTelomereLength: Double? = null

                if (!tumorOnly)
                {
                    // next we do the telomere length calculations from these files
                    val germlineTelLengthApp = TelLengthApp()
                    germlineTelLengthApp.params.sampleId = params.commonParams.referenceSampleId
                    germlineTelLengthApp.params.sampleType = SampleType.ref
                    germlineTelLengthApp.params.outputFile = germlineTelLegnthTsvPath()
                    germlineTelLengthApp.params.telbamFile = germlineTelbamPath()
                    germlineTelLengthApp.params.duplicatePercent = params.germlineDuplicateProportion
                    germlineTelLengthApp.params.meanReadsPerKb = params.germlineMeanReadsPerKb!!
                    germlineTelLengthApp.params.gc50ReadsPerKb = params.germlineGc50ReadsPerKb
                    germlineTelomereLength = germlineTelLengthApp.calcTelomereLength()
                }

                if (!germlineOnly)
                {
                    val tumorTelLengthApp = TelLengthApp()
                    tumorTelLengthApp.params.sampleId = params.commonParams.tumorSampleId
                    tumorTelLengthApp.params.sampleType = SampleType.tumor
                    tumorTelLengthApp.params.outputFile = tumorTelLengthTsvPath()
                    tumorTelLengthApp.params.telbamFile = tumorTelbamPath()
                    tumorTelLengthApp.params.germlineTelomereLength = germlineTelomereLength
                    tumorTelLengthApp.params.purity = params.tumorPurity
                    tumorTelLengthApp.params.ploidy = params.tumorPloidy
                    tumorTelLengthApp.params.duplicatePercent = params.tumorDuplicateProportion
                    tumorTelLengthApp.params.meanReadsPerKb = params.tumorMeanReadsPerKb!!
                    tumorTelLengthApp.params.gc50ReadsPerKb = params.tumorGc50ReadsPerKb
                    tumorTelLengthApp.calcTelomereLength()
                }

                // lastly we do the telomeric break ends
                if (!germlineOnly && !tumorOnly)
                {
                    val breakEndApp = BreakEndApp()
                    breakEndApp.params.sampleId = params.commonParams.tumorSampleId!!
                    breakEndApp.params.tumorTelbamFile = tumorTelbamPath()
                    breakEndApp.params.germlineTelbamFile = germlineTelbamPath()
                    breakEndApp.params.outputFile = breakEndTsvPath()
                    breakEndApp.params.refGenomeVersionStr = params.commonParams.getRefGenomeVersionStr()
                    breakEndApp.findBreakEnds()
                }
            }
            catch (e: InterruptedException)
            {
                logger.warn("Teal run interrupted, exiting")
                return 1
            }
            catch (e: IllegalStateException)
            {
                logger.error(e)
                return 1
            }

            val finish = Instant.now()
            val seconds = Duration.between(start, finish).seconds
            logger.info("Teal run complete, time taken: {}m {}s", seconds / 60, seconds % 60)
            return 0
        }

        fun germlineTelbamPath() : String
        {
            return "${params.commonParams.outputDir}/${params.commonParams.referenceSampleId}.teal.telbam.bam"
        }

        fun tumorTelbamPath() : String
        {
            return "${params.commonParams.outputDir}/${params.commonParams.tumorSampleId}.teal.telbam.bam"
        }

        fun germlineTelLegnthTsvPath() : String
        {
            return "${params.commonParams.outputDir}/${params.commonParams.referenceSampleId}.teal.tellength.tsv"
        }

        fun tumorTelLengthTsvPath() : String
        {
            return "${params.commonParams.outputDir}/${params.commonParams.tumorSampleId}.teal.tellength.tsv"
        }

        fun germlineTelReadTsvPath() : String
        {
            return "${params.commonParams.outputDir}/${params.commonParams.referenceSampleId}.teal.telread.tsv.gz"
        }

        fun tumorTelReadTsvPath() : String
        {
            return "${params.commonParams.outputDir}/${params.commonParams.tumorSampleId}.teal.telread.tsv.gz"
        }

        fun breakEndTsvPath() : String
        {
            return "${params.commonParams.outputDir}/${params.commonParams.tumorSampleId}.teal.breakend.tsv.gz"
        }
    }

    class PipelineMode
    {
        @ParametersDelegate
        val pipelineParams = TealPipelineParams()

        @ParametersDelegate
        val loggingOptions = LoggingOptions()

        fun run(): Int
        {
            pipelineParams.validate()
            loggingOptions.setLogLevel()
            val standaloneMode = StandaloneMode()
            standaloneMode.params = readPipelineFiles()
            return standaloneMode.run()
        }

        // from the pipeline locations
        private fun readPipelineFiles(): TealParams
        {
            val tealParams = TealParams(commonParams = pipelineParams.commonParams)

            if (pipelineParams.purple != null && pipelineParams.commonParams.tumorSampleId != null)
            {
                val purityContext = PurityContextFile.read(pipelineParams.purple!!, pipelineParams.commonParams.tumorSampleId!!)
                tealParams.tumorPurity = purityContext.bestFit().purity()
                tealParams.tumorPloidy = purityContext.bestFit().ploidy()
            }

            if (pipelineParams.cobalt != null)
            {
                if (pipelineParams.commonParams.tumorSampleId != null)
                {
                    val tumorGCMedianFilename =
                        GCMedianReadCountFile.generateFilename(pipelineParams.cobalt!!, pipelineParams.commonParams.tumorSampleId!!)
                    val tumorGCMedianReadCount = GCMedianReadCountFile.read(true, tumorGCMedianFilename)
                    tealParams.tumorMeanReadsPerKb = tumorGCMedianReadCount.meanReadCount()
                    tealParams.tumorGc50ReadsPerKb = tumorGCMedianReadCount.medianReadCount(ImmutableGCBucket(50))
                }

                if (pipelineParams.commonParams.referenceSampleId != null)
                {
                    val referenceGCMedianFilename =
                        GCMedianReadCountFile.generateFilename(pipelineParams.cobalt!!, pipelineParams.commonParams.referenceSampleId!!)
                    val referenceGCMedianReadCount = GCMedianReadCountFile.read(true, referenceGCMedianFilename)

                    tealParams.germlineMeanReadsPerKb = referenceGCMedianReadCount.meanReadCount()
                    tealParams.germlineGc50ReadsPerKb = referenceGCMedianReadCount.medianReadCount(ImmutableGCBucket(50))
                }
            }

            if (pipelineParams.tumorWgsMetrics != null)
            {
                // we need to try to guess the file name
                val metrics = WGSMetricsFile.read(pipelineParams.tumorWgsMetrics!!)
                logger.info("Loaded tumor WGS metrics from {}", pipelineParams.tumorWgsMetrics)

                tealParams.tumorDuplicateProportion = metrics.pctExcDupe()
            }

            if (pipelineParams.referenceWgsMetrics != null)
            {
                // we need to try to guess the file name
                val metrics = WGSMetricsFile.read(pipelineParams.referenceWgsMetrics!!)
                logger.info("Loaded reference WGS metrics from {}", pipelineParams.referenceWgsMetrics)
                tealParams.germlineDuplicateProportion = metrics.pctExcDupe()
            }

            logger.info("loaded teal params {} from pipeline files", tealParams)
            return tealParams
        }
    }

    companion object
    {
        private val logger = LogManager.getLogger(TealApplication::class.java)

        @JvmStatic
        fun main(args: Array<String>)
        {
            // here we have some voodoo to work out if we are being used in pipeline mode or the standalone mode
            val tealApp = StandaloneMode()
            val commanderStandalone = JCommander.newBuilder()
                .addObject(tealApp)
                .build()

            // use unix style formatter
            commanderStandalone.usageFormatter = UnixStyleUsageFormatter(commanderStandalone)
            commanderStandalone.parameterDescriptionComparator = DeclaredOrderParameterComparator(tealApp.javaClass)

            try
            {
                commanderStandalone.parse(*args)
                exitProcess(tealApp.run())
            }
            catch (standaloneParamException: ParameterException)
            {
                val tealPipelineApp = PipelineMode()
                val commanderPipeline = JCommander.newBuilder()
                    .addObject(tealPipelineApp)
                    .build()

                // use unix style formatter
                commanderPipeline.usageFormatter = UnixStyleUsageFormatter(commanderPipeline)
                commanderPipeline.parameterDescriptionComparator = DeclaredOrderParameterComparator(tealPipelineApp.javaClass)

                try
                {
                    commanderPipeline.parse(*args)
                    exitProcess(tealPipelineApp.run())
                }
                catch (pipelineParamException: ParameterException)
                {
                    println("Standalone mode:")
                    println("${standaloneParamException.message}")
                    commanderStandalone.usage()

                    println("Pipeline mode:")
                    println("${pipelineParamException.message}")
                    commanderPipeline.usage()

                    exitProcess(1)
                }
            }
        }
    }
}