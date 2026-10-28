package com.hartwig.hmftools.sigs.common;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.utils.MatrixFile;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class CommonUtils
{
    public static final String APP_NAME = "Sigs";

    public static final String SAMPLE_COUNTS_FILE = "sample_counts_file";
    public static final String SIGNATURES_FILE = "signatures_file";

    public static final String OUTPUT_FILE_ID = "output_file_id";

    public static final Logger SIG_LOGGER = LogManager.getLogger(CommonUtils.class);

    public static String parseOutputDir(final CommandLine cmd)
    {
        String outputDir = cmd.getOptionValue(OUTPUT_DIR);
        if(outputDir == null)
            return null;

        return checkAddDirSeparator(outputDir);
    }

    public static String formOutputFilename(final String outputDir, final String outputId, final String fileId)
    {
        if(outputId != null)
            return outputDir + fileId + "." + outputId + ".csv";
        else
            return outputDir + fileId  + ".csv";
    }

    public static Matrix loadSampleMatrixCounts(final String filename, final List<String> sampleIds)
    {
        return MatrixFile.loadMatrixDataFile(filename, sampleIds, Lists.newArrayList("BucketName"), false);
    }
}
