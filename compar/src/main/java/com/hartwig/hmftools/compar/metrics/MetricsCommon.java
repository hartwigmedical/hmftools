package com.hartwig.hmftools.compar.metrics;

public class MetricsCommon
{
    protected static final String FLD_MAPPED_PROPORTION = "MappedProportion";
    protected static final String FLD_DUPLICATE_PERCENTAGE = "DuplicatePercentage";

    protected static final double MAPPED_PROPORTION_ABS_THRESHOLD = 0.01;
    protected static final double MAPPED_PROPORTION_PCT_THRESHOLD = 0;
    protected static final double DUPLICATE_PERCENTAGE_ABS_THRESHOLD = 0.05;
    protected static final double DUPLICATE_PERCENTAGE_PCT_THRESHOLD = 0;

    protected static final String OLD_FLAGSTAT_FILE_EXTENSION = "_dedup.realigned.flagstat";
    protected static final String OLD_BAM_METRICS_FILE_EXTENSION = "_dedup_WGSMetrics.txt";
}
