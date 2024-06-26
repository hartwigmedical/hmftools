package com.hartwig.hmftools.teal.telbam

import com.beust.jcommander.JCommander
import com.beust.jcommander.ParametersDelegate
import com.beust.jcommander.UnixStyleUsageFormatter
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import com.hartwig.hmftools.common.utils.config.LoggingOptions
import com.hartwig.hmftools.common.utils.version.VersionInfo
import org.apache.logging.log4j.LogManager
import java.time.Duration
import java.time.Instant
import kotlin.system.exitProcess

class TelbamApp
{
    private val logger = LogManager.getLogger(javaClass)

    // add the BreakEndOptions to this command
    @ParametersDelegate
    val params = TelbamParams()

    @ParametersDelegate
    val loggingOptions = LoggingOptions()

    fun run(): Int
    {
        loggingOptions.setLogLevel()

        if (params.specificChromosomes.isNotEmpty())
        {
            logger.info("filtering for specific chromosomes: {}", params.specificChromosomes)
        }

        val versionInfo = VersionInfo("teal.version")
        logger.info("Teal version: {}", versionInfo.version())
        logger.info("starting telomeric analysis")
        val start = Instant.now()
        try
        {
            processBam()
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

    fun processBam()
    {
        logger.info("params: {}", params)
        val bamProcessor = BamProcessor(params)
        bamProcessor.processBam()
    }

    companion object
    {
        @JvmStatic
        fun main(args: Array<String>)
        {
            val telbamApp = TelbamApp()
            val commander = JCommander.newBuilder()
                .addObject(telbamApp)
                .build()

            // use unix style formatter
            commander.usageFormatter = UnixStyleUsageFormatter(commander)
            commander.parameterDescriptionComparator = DeclaredOrderParameterComparator(telbamApp.javaClass)

            try
            {
                commander.parse(*args)
            }
            catch (e: com.beust.jcommander.ParameterException)
            {
                println("${e.message}")
                commander.usage()
                exitProcess(1)
            }

            exitProcess(telbamApp.run())
        }
    }
}