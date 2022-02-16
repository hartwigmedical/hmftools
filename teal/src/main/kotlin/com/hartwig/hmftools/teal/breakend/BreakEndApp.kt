package com.hartwig.hmftools.teal.breakend

import com.beust.jcommander.JCommander
import com.hartwig.hmftools.common.utils.config.LoggingOptions
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import com.hartwig.hmftools.common.utils.version.VersionInfo
import org.apache.logging.log4j.LogManager
import java.time.Duration
import java.time.Instant
import kotlin.system.exitProcess

import com.beust.jcommander.ParametersDelegate
import com.beust.jcommander.UnixStyleUsageFormatter

class BreakEndApp
{
    private val logger = LogManager.getLogger(javaClass)

    // add the BreakEndOptions to this command
    @ParametersDelegate
    val params = BreakEndParams()

    @ParametersDelegate
    val loggingOptions = LoggingOptions()

    fun run(): Int
    {
        loggingOptions.setLogLevel()

        val versionInfo = VersionInfo("teal.version")
        logger.info("Teal version: {}", versionInfo.version())
        logger.info("starting telomeric breakend analysis")
        val start = Instant.now()

        findBreakEnds()

        val finish = Instant.now()
        val seconds = Duration.between(start, finish).seconds
        logger.info("Teal breakend run complete, time taken: {}m {}s", seconds / 60, seconds % 60)
        return 0
    }

    fun findBreakEnds()
    {
        val sampleBreakEndFinder = SampleBreakEndFinder(params)
        sampleBreakEndFinder.run()
    }

    companion object
    {
        @JvmStatic
        fun main(args: Array<String>)
        {
            val breakEndApp = BreakEndApp()
            val commander = JCommander.newBuilder()
                .addObject(breakEndApp)
                .build()

            // use unix style formatter
            commander.usageFormatter = UnixStyleUsageFormatter(commander)
            commander.parameterDescriptionComparator = DeclaredOrderParameterComparator(breakEndApp.javaClass)

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

            exitProcess(breakEndApp.run())
        }
    }
}
