package com.hartwig.hmftools.teal

import com.beust.jcommander.*
import com.hartwig.hmftools.common.genome.gc.GCMedianReadDepthFile
import com.hartwig.hmftools.common.genome.gc.ImmutableGCBucket
import com.hartwig.hmftools.common.metrics.BamMetricsSummary
import com.hartwig.hmftools.common.purple.PurityContextFile
import com.hartwig.hmftools.common.utils.config.LoggingOptions
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import org.apache.logging.log4j.LogManager
import kotlin.system.exitProcess

class TealPipelineApp
{
    @ParametersDelegate
    val pipelineParams = TealPipelineParams()

    @ParametersDelegate
    val loggingOptions = LoggingOptions()

    fun run(): Int
    {
        pipelineParams.validate()
        loggingOptions.setLogLevel()
        val standaloneMode = TealApplication()
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
                    GCMedianReadDepthFile.generateFilename(pipelineParams.cobalt!!, pipelineParams.commonParams.tumorSampleId!!)
                val tumorGCMedianReadDepth = GCMedianReadDepthFile.read(tumorGCMedianFilename)
                tealParams.tumorMeanReadDepth = tumorGCMedianReadDepth.meanReadDepth()
                tealParams.tumorGc50ReadDepth = tumorGCMedianReadDepth.medianReadDepth(ImmutableGCBucket(50))
            }

            if (pipelineParams.commonParams.referenceSampleId != null)
            {
                val referenceGCMedianFilename =
                    GCMedianReadDepthFile.generateFilename(pipelineParams.cobalt!!, pipelineParams.commonParams.referenceSampleId!!)
                val referenceGCMedianReadDepth = GCMedianReadDepthFile.read(referenceGCMedianFilename)

                tealParams.germlineMeanReadDepth = referenceGCMedianReadDepth.meanReadDepth()
                tealParams.germlineGc50ReadDepth = referenceGCMedianReadDepth.medianReadDepth(ImmutableGCBucket(50))
            }
        }

        if (pipelineParams.tumorWgsMetrics != null)
        {
            // we need to try to guess the file name
            val metrics = BamMetricsSummary.read(pipelineParams.tumorWgsMetrics!!)
            logger.info("Loaded tumor WGS metrics from {}", pipelineParams.tumorWgsMetrics)
            tealParams.tumorDuplicateProportion = metrics.duplicatePercent()
        }

        if (pipelineParams.referenceWgsMetrics != null)
        {
            // we need to try to guess the file name
            val metrics = BamMetricsSummary.read(pipelineParams.referenceWgsMetrics!!)
            logger.info("Loaded reference WGS metrics from {}", pipelineParams.referenceWgsMetrics)
            tealParams.germlineDuplicateProportion = metrics.duplicatePercent()
        }

        logger.info("loaded teal params {} from pipeline files", tealParams)
        return tealParams
    }

    companion object
    {
        private val logger = LogManager.getLogger(TealPipelineApp::class.java)

        @JvmStatic
        fun main(args: Array<String>)
        {
            // here we have some voodoo to work out if we are being used in pipeline mode or the standalone mode
            val tealApp = TealPipelineApp()
            val commander = JCommander.newBuilder()
                .addObject(tealApp)
                .build()

            // use unix style formatter
            commander.usageFormatter = UnixStyleUsageFormatter(commander)
            commander.parameterDescriptionComparator = DeclaredOrderParameterComparator(tealApp.javaClass)

            try
            {
                commander.parse(*args)
                exitProcess(tealApp.run())
            }
            catch (paramException: ParameterException)
            {
                println("${paramException.message}")
                commander.usage()
                exitProcess(1)
            }
        }
    }
}