package com.hartwig.hmftools.qsee.cohort;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

public class TrainConfig
{
    public final CommonPrepConfig CommonPrep;

    public final int NumPercentiles;
    public final double PercentileIncrement;

    private static final String NUM_PERCENTILES = "num_percentiles";
    private static final String NUM_PERCENTILES_DESC = "Number of percentiles (e.g. 5 would give percentiles: 0, 25, 50, 75, 100)";

    private static final String PERCENTILE_INCREMENT = "percentile_increment";
    private static final String PERCENTILE_INCREMENT_DESC = "Percentile increment (e.g. 25 would give percentiles: 0, 25, 50, 75, 100)";

    public TrainConfig(final ConfigBuilder configBuilder)
    {
        CommonPrep = new CommonPrepConfig(configBuilder);

        NumPercentiles = configBuilder.getInteger(NUM_PERCENTILES);
        PercentileIncrement = configBuilder.getDecimal(PERCENTILE_INCREMENT);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        CommonPrepConfig.registerConfig(configBuilder);

        configBuilder.addInteger(NUM_PERCENTILES, NUM_PERCENTILES_DESC, 101);
        configBuilder.addDecimal(PERCENTILE_INCREMENT, PERCENTILE_INCREMENT_DESC, Double.NaN);
    }

    public boolean hasPercentileIncrement() { return !Double.isNaN(PercentileIncrement); }
}

