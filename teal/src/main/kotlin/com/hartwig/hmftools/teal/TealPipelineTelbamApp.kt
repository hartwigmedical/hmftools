package com.hartwig.hmftools.teal

import com.beust.jcommander.*
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import com.hartwig.hmftools.common.utils.config.LoggingOptions
import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import com.hartwig.hmftools.common.utils.version.VersionInfo
import com.hartwig.hmftools.teal.telbam.TelbamApp
import org.apache.logging.log4j.LogManager
import java.time.Duration
import java.time.Instant
import kotlin.system.exitProcess

// this is special mode which only process the bam and create telbams for the
// pipeline run. The reason for splitting into two part is for performance reasons.
// It allows us to avoid reprocessing the full bam.
class TealPipelineTelbamApp
{
    @ParametersDelegate
    val tealCommonParams = TealCommonParams()

    @ParametersDelegate
    val loggingOptions = LoggingOptions()

    fun run(): Int
    {
        tealCommonParams.validate()
        loggingOptions.setLogLevel()

        if (!FileWriterUtils.checkCreateOutputDir(tealCommonParams.outputDir))
        {
            logger.error("failed to create output directory({})", tealCommonParams.outputDir)
            return 1
        }

        val versionInfo = VersionInfo("teal.version")
        logger.info("Teal version: {}", versionInfo.version())
        logger.info("starting telomeric analysis")
        logger.info("{}", tealCommonParams)
        val start = Instant.now()

        val tumorOnly = tealCommonParams.tumorOnly()
        val germlineOnly = tealCommonParams.referenceOnly()

        try
        {
            logger.info("creating telbam files")

            // first we generate the telbam files
            // we pass the params to the telbam app

            if (!tumorOnly)
            {
                val germlineTelbamApp = TelbamApp()
                germlineTelbamApp.params.bamFile = tealCommonParams.referenceBamFile!!
                germlineTelbamApp.params.refGenomeFile = tealCommonParams.refGenomeFile
                germlineTelbamApp.params.telbamFile = tealCommonParams.germlineTelbamPath()
                germlineTelbamApp.params.threadCount = tealCommonParams.threadCount
                germlineTelbamApp.processBam()
            }

            if (!germlineOnly)
            {
                val tumorTelbamApp = TelbamApp()
                tumorTelbamApp.params.bamFile = tealCommonParams.tumorBamFile!!
                tumorTelbamApp.params.refGenomeFile = tealCommonParams.refGenomeFile
                tumorTelbamApp.params.telbamFile = tealCommonParams.tumorTelbamPath()
                tumorTelbamApp.params.threadCount = tealCommonParams.threadCount
                tumorTelbamApp.processBam()
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
    
    companion object
    {
        private val logger = LogManager.getLogger(TealPipelineTelbamApp::class.java)

        @JvmStatic
        fun main(args: Array<String>)
        {
            // here we have some voodoo to work out if we are being used in pipeline mode or the standalone mode
            val app = TealPipelineTelbamApp()
            val jCommander = JCommander.newBuilder()
                .addObject(app)
                .build()

            // use unix style formatter
            jCommander.usageFormatter = UnixStyleUsageFormatter(jCommander)
            jCommander.parameterDescriptionComparator = DeclaredOrderParameterComparator(app.javaClass)

            try
            {
                jCommander.parse(*args)
                exitProcess(app.run())
            }
            catch (paramException: ParameterException)
            {
                println("${paramException.message}")
                jCommander.usage()
                exitProcess(1)
            }
        }
    }
}