package com.hartwig.hmftools.teal

import com.hartwig.hmftools.common.utils.config.ConfigBuilder
import com.hartwig.hmftools.common.utils.config.ConfigUtils
import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import com.hartwig.hmftools.common.utils.version.VersionInfo
import com.hartwig.hmftools.teal.breakend.BreakEndApp
import com.hartwig.hmftools.teal.telbam.TelbamApp
import com.hartwig.hmftools.teal.tellength.SampleType
import com.hartwig.hmftools.teal.tellength.TelLengthApp
import org.apache.logging.log4j.LogManager
import java.time.Duration
import java.time.Instant
import kotlin.system.exitProcess

class TealApplication(val params: TealParams)
{
    fun run(): Int
    {
        params.validate()

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
                germlineTelbamApp.params.telbamFile = params.commonParams.germlineTelbamPath()
                germlineTelbamApp.params.threadCount = params.commonParams.threadCount
                germlineTelbamApp.processBam()
            }

            if (!germlineOnly)
            {
                val tumorTelbamApp = TelbamApp()
                tumorTelbamApp.params.bamFile = params.commonParams.tumorBamFile!!
                tumorTelbamApp.params.refGenomeFile = params.commonParams.refGenomeFile
                tumorTelbamApp.params.telbamFile = params.commonParams.tumorTelbamPath()
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
                germlineTelLengthApp.params.outputFile = params.commonParams.germlineTelLegnthTsvPath()
                germlineTelLengthApp.params.telbamFile = params.commonParams.germlineTelbamPath()
                germlineTelLengthApp.params.duplicatePercent = params.germlineDuplicateProportion
                germlineTelLengthApp.params.meanReadDepth = params.germlineMeanReadDepth!!
                germlineTelLengthApp.params.gc50ReadDepth = params.germlineGc50ReadDepth
                germlineTelomereLength = germlineTelLengthApp.calcTelomereLength()
            }

            if (!germlineOnly)
            {
                val tumorTelLengthApp = TelLengthApp()
                tumorTelLengthApp.params.sampleId = params.commonParams.tumorSampleId
                tumorTelLengthApp.params.sampleType = SampleType.tumor
                tumorTelLengthApp.params.outputFile = params.commonParams.tumorTelLengthTsvPath()
                tumorTelLengthApp.params.telbamFile = params.commonParams.tumorTelbamPath()
                tumorTelLengthApp.params.germlineTelomereLength = germlineTelomereLength
                tumorTelLengthApp.params.purity = params.tumorPurity
                tumorTelLengthApp.params.ploidy = params.tumorPloidy
                tumorTelLengthApp.params.duplicatePercent = params.tumorDuplicateProportion
                tumorTelLengthApp.params.meanReadDepth = params.tumorMeanReadDepth!!
                tumorTelLengthApp.params.gc50ReadDepth = params.tumorGc50ReadDepth
                tumorTelLengthApp.calcTelomereLength()
            }

            // lastly we do the telomeric break ends
            if (!germlineOnly && !tumorOnly)
            {
                val breakEndApp = BreakEndApp()
                breakEndApp.params.sampleId = params.commonParams.tumorSampleId!!
                breakEndApp.params.tumorTelbamFile = params.commonParams.tumorTelbamPath()
                breakEndApp.params.germlineTelbamFile = params.commonParams.germlineTelbamPath()
                breakEndApp.params.outputFile = params.commonParams.breakEndTsvPath()
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

    companion object
    {
        private val logger = LogManager.getLogger(TealApplication::class.java)

        @JvmStatic
        fun main(args: Array<String>)
        {
            val configBuilder = ConfigBuilder("Teal")
            TealParams.registerConfig(configBuilder)
            ConfigUtils.addLoggingOptions(configBuilder)
            configBuilder.checkAndParseCommandLine(args)
            val tealApp = TealApplication(TealParams(configBuilder))
            exitProcess(tealApp.run())
        }
    }
}