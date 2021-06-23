package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_MAX_ALT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_MIN_MAP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_SAMPLE_SIZE;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class QualityRecalibrationConfig
{
    public final boolean Enabled;
    public final boolean Plot;
    public final int MaxAltCount;
    public final int SampleSize;
    public final int MinMapQuality;

    private static final String BQR_PLOT = "bqr_plot";
    private static final String BQR_ENABLED = "bqr_enabled";
    private static final String BQR_SAMPLE_SIZE = "bqr_sample_size";
    private static final String BQR_MAX_ALT_COUNT = "bqr_max_alt_count";
    private static final String BQR_MIN_MAP_QUAL = "bqr_min_map_qual";

    private static final boolean DEFAULT_BQR_PLOT = true;
    private static final boolean DEFAULT_BQR_ENABLED = true;

    public QualityRecalibrationConfig(final CommandLine cmd)
    {
        Enabled = getConfigValue(cmd, BQR_ENABLED, DEFAULT_BQR_ENABLED);
        Plot = getConfigValue(cmd, BQR_PLOT, DEFAULT_BQR_PLOT);
        MaxAltCount = getConfigValue(cmd, BQR_MAX_ALT_COUNT, DEFAULT_BQR_MAX_ALT_COUNT);
        SampleSize = getConfigValue(cmd, BQR_SAMPLE_SIZE, DEFAULT_BQR_SAMPLE_SIZE);
        MinMapQuality = getConfigValue(cmd, BQR_MIN_MAP_QUAL, DEFAULT_BQR_MIN_MAP_QUAL);
    }

    public QualityRecalibrationConfig()
    {
        Enabled = false;
        Plot = false;
        MaxAltCount = DEFAULT_BQR_MAX_ALT_COUNT;
        SampleSize = DEFAULT_BQR_SAMPLE_SIZE;
        MinMapQuality = DEFAULT_BQR_MIN_MAP_QUAL;
    }

    @NotNull
    public static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(BQR_ENABLED, true, "BQR (Base Quality Recalibration) enabled [" + DEFAULT_BQR_ENABLED + "]");
        options.addOption(BQR_PLOT, true, "BQR plots [" + DEFAULT_BQR_PLOT + "]");
        options.addOption(BQR_MAX_ALT_COUNT, true, "BQR maximum alt count to be an error [" + DEFAULT_BQR_MAX_ALT_COUNT + "]");
        options.addOption(BQR_SAMPLE_SIZE, true, "BQR sampling size per autosome [" + DEFAULT_BQR_SAMPLE_SIZE + "]");
        options.addOption(BQR_MIN_MAP_QUAL, true, "BQR min base quality remap qual [" + DEFAULT_BQR_MIN_MAP_QUAL + "]");
        return options;
    }
}
