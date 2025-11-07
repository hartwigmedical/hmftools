package com.hartwig.hmftools.qsee.vis;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

public class VisPrepConfig
{
    public final CommonPrepConfig CommonPrep;

    public final String CohortPercentilesFile;

    public static final String COHORT_PERCENTILES_FILE = "cohort_percentiles_file";
    public static final String COHORT_PERCENTILES_FILE_DESC = "Path to the cohort percentiles file";

    public VisPrepConfig(final ConfigBuilder configBuilder)
    {
        CommonPrep = new CommonPrepConfig(configBuilder);

        CohortPercentilesFile = configBuilder.getValue(COHORT_PERCENTILES_FILE);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        CommonPrepConfig.registerConfig(configBuilder);

        configBuilder.addPath(COHORT_PERCENTILES_FILE, false, COHORT_PERCENTILES_FILE_DESC);
    }
}
