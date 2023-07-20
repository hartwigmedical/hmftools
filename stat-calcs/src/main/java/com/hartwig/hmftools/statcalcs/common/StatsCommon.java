package com.hartwig.hmftools.statcalcs.common;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class StatsCommon
{
    public static final String OUTPUT_FILE_ID = "output_file_id";

    public static final Logger STAT_LOGGER = LogManager.getLogger(StatsCommon.class);

    public static String formOutputFilename(final String outputDir, final String outputId, final String fileId)
    {
        if(outputId != null)
            return outputDir + fileId + "." + outputId + ".csv";
        else
            return outputDir + fileId  + ".csv";
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }

}
