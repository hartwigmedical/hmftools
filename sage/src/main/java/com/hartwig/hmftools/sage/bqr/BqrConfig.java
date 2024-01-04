package com.hartwig.hmftools.sage.bqr;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_MAX_ALT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_MAX_ALT_PERC;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_MIN_MAP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_SAMPLE_SIZE;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class BqrConfig
{
    public final boolean Enabled;
    public final boolean LoadBqrFiles;
    public final boolean WriteFile;
    public final boolean WritePlot;
    public final double MaxAltPerc;
    public final int MaxAltCount;
    public final int SampleSize;
    public final int MinMapQuality;

    private static final String DISABLE_BQR = "disable_bqr";
    private static final String LOAD_BQR_FILES = "load_bqr";
    private static final String WRITE_BQR_DATA = "write_bqr_data";
    private static final String WRITE_BQR_PLOT = "write_bqr_plot";

    private static final String BQR_SAMPLE_SIZE = "bqr_sample_size";
    private static final String BQR_MAX_ALT_PERC = "bqr_max_alt_perc";
    private static final String BQR_MAX_ALT_COUNT = "bqr_max_alt_count";
    private static final String BQR_MIN_MAP_QUAL = "bqr_min_map_qual";

    public BqrConfig(final ConfigBuilder configBuilder)
    {
        Enabled = !configBuilder.hasFlag(DISABLE_BQR);

        if(configBuilder.hasFlag(LOAD_BQR_FILES))
        {
            LoadBqrFiles = true;
            WriteFile = false;
            WritePlot = false;
        }
        else
        {
            LoadBqrFiles = false;
            WriteFile = Enabled || configBuilder.hasFlag(WRITE_BQR_DATA); // written by default now
            WritePlot = configBuilder.hasFlag(WRITE_BQR_PLOT);
        }

        MaxAltPerc = configBuilder.getDecimal(BQR_MAX_ALT_PERC);
        MaxAltCount = configBuilder.getInteger(BQR_MAX_ALT_COUNT);
        SampleSize = configBuilder.getInteger(BQR_SAMPLE_SIZE);
        MinMapQuality = configBuilder.getInteger(BQR_MIN_MAP_QUAL);
    }

    public BqrConfig()
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

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(DISABLE_BQR, "Disable Base Quality Recalibration");
        configBuilder.addFlag(WRITE_BQR_DATA, "Write BQR output file");
        configBuilder.addFlag(WRITE_BQR_PLOT, "Generate BQR plot");
        configBuilder.addFlag(LOAD_BQR_FILES, "Attemps to find and load previously-written BQR files");
        configBuilder.addDecimal(BQR_MAX_ALT_PERC, "BQR maximum alt percent to be an error", DEFAULT_BQR_MAX_ALT_PERC);
        configBuilder.addInteger(BQR_MAX_ALT_COUNT, "BQR maximum alt count to be an error", DEFAULT_BQR_MAX_ALT_COUNT);
        configBuilder.addInteger(BQR_SAMPLE_SIZE, "BQR sampling size per autosome", DEFAULT_BQR_SAMPLE_SIZE);
        configBuilder.addInteger(BQR_MIN_MAP_QUAL, "BQR min base quality remap qual", DEFAULT_BQR_MIN_MAP_QUAL);
    }
}
