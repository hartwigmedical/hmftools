package com.hartwig.hmftools.teal.breakend

import com.hartwig.hmftools.common.utils.FileWriterUtils
import com.hartwig.hmftools.common.utils.ConfigUtils
import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import org.apache.logging.log4j.LogManager

const val DEFAULT_BREAK_POINT_MERGE_DISTANCE = 50

class BreakEndAppConfig(cmd: CommandLine)
{
    val SampleId: String
    val TumorTelbamFile: String
    val GermlineTelbamFile: String
    val OutputDir: String
    val SplitTelomereMatchThreshold: Double = 0.9

    // I find using 0.8 can still get all the break ends gripss found.
    // using 0.7 and we start excluding some
    val AlignedTelomereMatchThreshold: Double = 0.8

    val BreakPointMergeDistance: Int

    fun isValid(): Boolean
    {
        if (!FileWriterUtils.checkCreateOutputDir(OutputDir))
        {
            LOGGER.error("failed to create output directory({})", OutputDir)
            return false
        }
        return true
    }

    companion object
    {
        private const val SAMPLE_ID = "sample_id"
        private const val TUMOR_TELBAM_FILE = "tumor_telbam_file"
        private const val GERMLINE_TELBAM_FILE = "germline_telbam_file"
        private const val BOUNDARY_ZONE = "boundary_zone"
        private const val BREAKPOINT_MERGE_DISTANCE = "breakpoint_merge_distance"
        val LOGGER = LogManager.getLogger(BreakEndAppConfig::class.java)
        fun createOptions(): Options
        {
            val options = Options()
            options.addOption(Option.builder(SAMPLE_ID).hasArg().required().desc("ID of tumor sample").build())
            options.addOption(Option.builder(TUMOR_TELBAM_FILE).hasArg().required().desc("Path to tumor telbam.bam file").build())
            options.addOption(Option.builder(GERMLINE_TELBAM_FILE).hasArg().required().desc("Path to reference telbam.bam file").build())
            options.addOption(Option.builder(FileWriterUtils.OUTPUT_DIR).hasArg().required().desc("Output directory").build())
            options.addOption(Option.builder(BREAKPOINT_MERGE_DISTANCE).hasArg().type(Number::class.java).
                    desc("breakpoints within this distance will be merged (default = $DEFAULT_BREAK_POINT_MERGE_DISTANCE").build())
            ConfigUtils.addLoggingOptions(options)
            return options
        }
    }

    init
    {
        SampleId = cmd.getOptionValue(SAMPLE_ID)
        TumorTelbamFile = cmd.getOptionValue(TUMOR_TELBAM_FILE)
        GermlineTelbamFile = cmd.getOptionValue(GERMLINE_TELBAM_FILE)
        OutputDir = FileWriterUtils.parseOutputDir(cmd)
        BreakPointMergeDistance = cmd.getOptionValue(BREAKPOINT_MERGE_DISTANCE, DEFAULT_BREAK_POINT_MERGE_DISTANCE.toString()).toInt()
    }
}