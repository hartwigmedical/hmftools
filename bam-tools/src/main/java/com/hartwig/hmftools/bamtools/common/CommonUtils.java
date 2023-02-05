package com.hartwig.hmftools.bamtools.common;

import static com.hartwig.hmftools.bamtools.metrics.MetricsConfig.BT_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;

public final class CommonUtils
{
    // config strings
    public static final String SAMPLE = "sample";
    public static final String BAM_FILE = "bam_file";
    public static final String PARTITION_SIZE = "partition_size";

    public static final int DEFAULT_CHR_PARTITION_SIZE = 1000000;

    public static boolean checkFileExists(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
        {
            BT_LOGGER.error("invalid file path: {}", filename);
            return false;
        }

        return true;
    }

    public static String formFilename(final String sampleId, final String outputDir, final String outputId, final String fileType)
    {
        String filename = outputDir + sampleId;

        filename += ".bam_" + fileType;

        if(outputId != null)
            filename += "." + outputId;

        filename += ".csv";

        return filename;
    }
}
