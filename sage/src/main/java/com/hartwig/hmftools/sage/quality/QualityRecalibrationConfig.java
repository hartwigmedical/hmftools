package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_MAX_ALT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_MAX_ALT_PERC;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_MIN_MAP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_SAMPLE_SIZE;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class QualityRecalibrationConfig
{
    public final boolean Enabled;
    public final boolean LoadBqrFiles;
    public final boolean WriteFile;
    public final boolean WritePlot;
    public final double MaxAltPerc;
    public final int MaxAltCount;
    public final int SampleSize;
    public final int MinMapQuality;

    private static final String BQR_ENABLED = "bqr_enabled";
    private static final String BQR_SAMPLE_SIZE = "bqr_sample_size";
    private static final String BQR_MAX_ALT_PERC = "bqr_max_alt_perc";
    private static final String BQR_MAX_ALT_COUNT = "bqr_max_alt_count";
    private static final String BQR_MIN_MAP_QUAL = "bqr_min_map_qual";
    private static final String WRITE_BQR_DATA = "write_bqr_data";
    private static final String WRITE_BQR_PLOT = "write_bqr_plot";
    private static final String LOAD_BQR_FILES = "load_bqr_files";

    private static final boolean DEFAULT_BQR_ENABLED = true;

    public QualityRecalibrationConfig(final CommandLine cmd)
    {
        Enabled = getConfigValue(cmd, BQR_ENABLED, DEFAULT_BQR_ENABLED);

        if(cmd.hasOption(LOAD_BQR_FILES))
        {
            LoadBqrFiles = true;
            WriteFile = false;
            WritePlot = false;
        }
        else
        {
            LoadBqrFiles = false;
            WriteFile = cmd.hasOption(WRITE_BQR_DATA);
            WritePlot = cmd.hasOption(WRITE_BQR_PLOT);
        }

        MaxAltPerc = getConfigValue(cmd, BQR_MAX_ALT_PERC, DEFAULT_BQR_MAX_ALT_PERC);
        MaxAltCount = getConfigValue(cmd, BQR_MAX_ALT_COUNT, DEFAULT_BQR_MAX_ALT_COUNT);
        SampleSize = getConfigValue(cmd, BQR_SAMPLE_SIZE, DEFAULT_BQR_SAMPLE_SIZE);
        MinMapQuality = getConfigValue(cmd, BQR_MIN_MAP_QUAL, DEFAULT_BQR_MIN_MAP_QUAL);
    }

    public QualityRecalibrationConfig()
    {
        Enabled = false;
        WritePlot = false;
        LoadBqrFiles = false;
        WriteFile = false;
        MaxAltPerc = DEFAULT_BQR_MAX_ALT_PERC;
        MaxAltCount = DEFAULT_BQR_MAX_ALT_COUNT;
        SampleSize = DEFAULT_BQR_SAMPLE_SIZE;
        MinMapQuality = DEFAULT_BQR_MIN_MAP_QUAL;
    }

    @NotNull
    public static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(BQR_ENABLED, true, "BQR (Base Quality Recalibration) enabled [" + DEFAULT_BQR_ENABLED + "]");
        options.addOption(WRITE_BQR_DATA, false, "Write BQR output file");
        options.addOption(WRITE_BQR_PLOT, false, "Generate BQR plot");
        options.addOption(LOAD_BQR_FILES, false, "Attemps to find and load previously-written BQR files");
        options.addOption(BQR_MAX_ALT_PERC, true, "BQR maximum alt percent to be an error [" + DEFAULT_BQR_MAX_ALT_PERC + "]");
        options.addOption(BQR_MAX_ALT_COUNT, true, "BQR maximum alt count to be an error [" + DEFAULT_BQR_MAX_ALT_COUNT + "]");
        options.addOption(BQR_SAMPLE_SIZE, true, "BQR sampling size per autosome [" + DEFAULT_BQR_SAMPLE_SIZE + "]");
        options.addOption(BQR_MIN_MAP_QUAL, true, "BQR min base quality remap qual [" + DEFAULT_BQR_MIN_MAP_QUAL + "]");
        return options;
    }
}
