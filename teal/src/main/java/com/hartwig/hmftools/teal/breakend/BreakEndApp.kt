package com.hartwig.hmftools.teal.breakend

import com.hartwig.hmftools.common.utils.ConfigUtils
import com.hartwig.hmftools.common.utils.version.VersionInfo
import com.hartwig.hmftools.teal.TeloConfig
import org.apache.commons.cli.*
import java.time.Duration
import java.time.Instant
import kotlin.system.exitProcess

class BreakEndApp(options: Options, args: Array<String>)
{
    private val mConfig: BreakEndAppConfig
    fun run()
    {
        if (!mConfig.isValid())
        {
            TeloConfig.TE_LOGGER.error(" invalid config, exiting")
            exitProcess(1)
        }
        TeloConfig.TE_LOGGER.info("starting telomeric analysis")
        val start = Instant.now()

        val sampleBreakEndFinder = SampleBreakEndFinder(mConfig)
        sampleBreakEndFinder.run()

        val finish = Instant.now()
        val seconds = Duration.between(start, finish).seconds
        TeloConfig.TE_LOGGER.info("Telo run complete, time taken: {}m {}s", seconds / 60, seconds % 60)
    }

    init
    {
        val versionInfo = VersionInfo("teal.version")
        TeloConfig.TE_LOGGER.info("Telo version: {}", versionInfo.version())
        val cmd = createCommandLine(args, options)
        ConfigUtils.setLogLevel(cmd)
        mConfig = BreakEndAppConfig(cmd)
    }
}

fun main(args: Array<String>)
{
    val options = BreakEndAppConfig.createOptions()
    val application : BreakEndApp
    try
    {
        application = BreakEndApp(options, args)
        application.run()
    } catch (e: ParseException)
    {
        TeloConfig.TE_LOGGER.warn(e)
        val formatter = HelpFormatter()
        formatter.printHelp("TeloApplication", options)
        exitProcess(1)
    }
}

@Throws(ParseException::class)
private fun createCommandLine(args: Array<String>, options: Options): CommandLine
{
    val parser: CommandLineParser = DefaultParser()
    return parser.parse(options, args)
}
