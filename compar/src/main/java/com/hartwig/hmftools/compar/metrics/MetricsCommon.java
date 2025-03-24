package com.hartwig.hmftools.compar.metrics;

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

    public static BamMetricsSummary loadBamMetricsSummary(final String sampleId, String directory) throws IOException
    {
        String currentFileName = BamMetricsSummary.generateFilename(directory, sampleId);
        return BamMetricsSummary.read(currentFileName);
    }

    public static String determineFlagStatsFilePath(final String sampleId, final String directory)
    {
        return BamFlagStats.generateFilename(directory, sampleId);
    }
}
