package com.hartwig.hmftools.teal.tellength

import com.beust.jcommander.JCommander
import com.beust.jcommander.ParametersDelegate
import com.beust.jcommander.UnixStyleUsageFormatter

import com.hartwig.hmftools.common.utils.config.LoggingOptions
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import com.hartwig.hmftools.common.utils.version.VersionInfo
import com.hartwig.hmftools.teal.telbam.TelbamReader
import org.apache.logging.log4j.LogManager
import java.io.File
import java.time.Duration
import java.time.Instant
import kotlin.system.exitProcess

class TelLengthApp
{
    private val logger = LogManager.getLogger(javaClass)

    @ParametersDelegate
    val params = TelLengthParams()

    @ParametersDelegate
    val loggingOptions = LoggingOptions()

    fun run(): Int
    {
        loggingOptions.setLogLevel()

        if (params.sampleId == null)
        {
            logger.error("sample Id cannot be null")
            return 1
        }

        if (params.sampleType == null)
        {
            logger.error("sample type cannot be null")
            return 1
        }

        if (params.telbamFile == null)
        {
            logger.error("telbam file path cannot be null")
            return 1
        }

        val versionInfo = VersionInfo("teal.version")
        logger.info("Teal version: {}", versionInfo.version())
        logger.info("starting telomere length analysis")

        val start = Instant.now()

        if (params.gc50ReadsPerKb == null)
        {
            logger.info("no gc50ReadsPerKb provided, using meanReadsPerKb({})", params.meanReadsPerKb)
            params.gc50ReadsPerKb = params.meanReadsPerKb
        }

        calcTelomereLength()

        val finish = Instant.now()
        val seconds = Duration.between(start, finish).seconds
        logger.info("Teal telomere length run complete, time taken: {}m {}s", seconds / 60, seconds % 60)
        return 0
    }

    fun calcTelomereLength() : Double
    {
        logger.info("params: {}", params)
        logger.info("processing telbam: {}", params.telbamFile)
        val telbamReader = TelbamReader(File(params.telbamFile!!))
        telbamReader.read()

        val analyser = TelomereLengthCalc(
            purity = params.purity,
            ploidy = params.ploidy,
            duplicateProportion = params.duplicatePercent,
            meanReadsPerKb = params.meanReadsPerKb,
            gc50ReadsPerKb = params.gc50ReadsPerKb!!,
            germlineTelomereLength = params.germlineTelomereLength)
        analyser.processReadGroups(telbamReader.readGroups.values)
        analyser.writeTelLengthTsv(params.outputFile!!, params.sampleId!!, params.sampleType!!)
        return analyser.calcTelomereLength()
    }

    companion object
    {
        @JvmStatic
        fun main(args: Array<String>)
        {
            val telLengthApp = TelLengthApp()
            val commander = JCommander.newBuilder()
                .addObject(telLengthApp)
                .build()

            // use unix style formatter
            commander.usageFormatter = UnixStyleUsageFormatter(commander)
            commander.parameterDescriptionComparator = DeclaredOrderParameterComparator(telLengthApp.javaClass)

            try
            {
                commander.parse(*args)
                telLengthApp.params.validate()
            }
            catch (e: com.beust.jcommander.ParameterException)
            {
                println("${e.message}")
                commander.usage()
                exitProcess(1)
            }

            exitProcess(telLengthApp.run())
        }
    }
}