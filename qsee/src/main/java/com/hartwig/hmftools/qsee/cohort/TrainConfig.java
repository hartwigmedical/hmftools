package com.hartwig.hmftools.qsee.cohort;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.prep.PrepConfig;

public class TrainConfig
{
    public final PrepConfig Prep;

    public final int NumPercentiles;
    public final double PercentileInterval;

    private static final String NUM_PERCENTILES = "num_percentiles";
    private static final String NUM_PERCENTILES_DESC = "Number of percentiles (e.g. 5 would give percentiles: 0, 25, 50, 75, 100)";

    private static final String PERCENTILE_INTERVAL = "percentile_interval";
    private static final String PERCENTILE_INTERVAL_DESC = "Percentile interval (e.g. 25 would give percentiles: 0, 25, 50, 75, 100)";

    public TrainConfig(final ConfigBuilder configBuilder)
    {
        Prep = new PrepConfig(configBuilder);
        NumPercentiles = configBuilder.getInteger(NUM_PERCENTILES);
        PercentileInterval = configBuilder.getDecimal(PERCENTILE_INTERVAL);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        PrepConfig.registerConfig(configBuilder);

        configBuilder.addInteger(NUM_PERCENTILES, NUM_PERCENTILES_DESC, 101);
        configBuilder.addDecimal(PERCENTILE_INTERVAL, PERCENTILE_INTERVAL_DESC, Double.NaN);
    }

    public boolean hasPercentileInterval() { return !Double.isNaN(PercentileInterval); }
}

