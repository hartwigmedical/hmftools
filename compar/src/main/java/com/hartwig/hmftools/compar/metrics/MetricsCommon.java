package com.hartwig.hmftools.compar.metrics;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.common.CommonUtils.fileExists;

import java.io.IOException;

import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricsSummary;

public class MetricsCommon
{
    protected static final String FLD_MAPPED_PROPORTION = "MappedProportion";
    protected static final String FLD_DUPLICATE_PERCENTAGE = "DuplicatePercentage";

    protected static final double MAPPED_PROPORTION_ABS_THRESHOLD = 0.01;
    protected static final double MAPPED_PROPORTION_PCT_THRESHOLD = 0;
    protected static final double DUPLICATE_PERCENTAGE_ABS_THRESHOLD = 0.05;
    protected static final double DUPLICATE_PERCENTAGE_PCT_THRESHOLD = 0;

    protected static final String OLD_FLAGSTAT_FILE_EXTENSION = ".flagstat";
    protected static final String VERY_OLD_FLAGSTAT_FILE_EXTENSION = "_dedup.realigned.flagstat";
    protected static final String OLD_BAM_METRICS_FILE_EXTENSION = ".wgsmetrics";
    protected static final String VERY_OLD_BAM_METRICS_FILE_EXTENSION = "_dedup_WGSMetrics.txt";

    public static BamMetricsSummary loadBamMetricsSummary(final String sampleId, String directory) throws IOException
    {
        String currentFileName = BamMetricsSummary.generateFilename(directory, sampleId);
        String oldFileName = checkAddDirSeparator(directory) + sampleId + OLD_BAM_METRICS_FILE_EXTENSION;
        String veryOldFileName = checkAddDirSeparator(directory) + sampleId + VERY_OLD_BAM_METRICS_FILE_EXTENSION;
        if(!fileExists(currentFileName) && fileExists(oldFileName))
        {
            return OldWGSMetricsFile.read(oldFileName);
        }
        else if(!fileExists(currentFileName) && !fileExists(oldFileName) && fileExists(veryOldFileName))
        {
            return OldWGSMetricsFile.read(veryOldFileName);
        }
        else
        {
            return BamMetricsSummary.read(currentFileName);
        }
    }

    public static String determineFlagStatsFilePath(final String sampleId, final String directory)
    {
        String currentFileName = BamFlagStats.generateFilename(directory, sampleId);
        String oldFileName = checkAddDirSeparator(directory) + sampleId + OLD_FLAGSTAT_FILE_EXTENSION;
        String veryOldFileName = checkAddDirSeparator(directory) + sampleId + VERY_OLD_FLAGSTAT_FILE_EXTENSION;
        if(!fileExists(currentFileName) && fileExists(oldFileName))
        {
            return oldFileName;
        }
        else if(!fileExists(currentFileName) && !fileExists(oldFileName) && fileExists(veryOldFileName))
        {
            return veryOldFileName;
        }
        else
        {
            return currentFileName;
        }
    }
}
