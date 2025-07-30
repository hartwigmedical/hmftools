package com.hartwig.hmftools.redux.bqr;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class BqrConfig
{
    public final boolean Enabled;
    public final boolean UseAllRegions;
    public final boolean WritePlot;

    private static final String DISABLE_BQR = "bqr_disable";
    private static final String WRITE_BQR_PLOT = "bqr_write_plot";
    private static final String USE_ALL_REGIONS = "bqr_use_all_regions";

    public BqrConfig(final ConfigBuilder configBuilder)
    {
        Enabled = !configBuilder.hasFlag(DISABLE_BQR);
        WritePlot = configBuilder.hasFlag(WRITE_BQR_PLOT);
        UseAllRegions = configBuilder.hasFlag(USE_ALL_REGIONS);
    }

    public BqrConfig()
    {
        Enabled = false;
        WritePlot = false;
        UseAllRegions = false;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(DISABLE_BQR, "Disable Base Quality Recalibration");
        configBuilder.addFlag(WRITE_BQR_PLOT, "Generate BQR plot");
        configBuilder.addFlag(USE_ALL_REGIONS, "Run BQR from full BAM");
    }
}
