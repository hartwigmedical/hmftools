package com.hartwig.hmftools.statcalcs.common;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;

import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class StatsCommon
{
    public static final String LOG_DEBUG = "log_debug";

    public static final String OUTPUT_FILE_ID = "output_file_id";

    public static final Logger STAT_LOGGER = LogManager.getLogger(StatsCommon.class);

    public static String formOutputFilename(final String outputDir, final String outputId, final String fileId)
    {
        if(outputId != null)
            return outputDir + fileId + "." + outputId + ".csv";
        else
            return outputDir + fileId  + ".csv";
    }

    public static void addCmdLineArgs(final Options options)
    {
        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }

}
