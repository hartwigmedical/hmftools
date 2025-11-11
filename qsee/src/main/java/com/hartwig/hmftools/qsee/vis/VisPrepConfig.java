package com.hartwig.hmftools.qsee.vis;

import static com.hartwig.hmftools.qsee.cohort.CohortPercentilesConfig.COHORT_PERCENTILES_FILE;
import static com.hartwig.hmftools.qsee.cohort.CohortPercentilesConfig.COHORT_PERCENTILES_FILE_DESC;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

public class VisPrepConfig
{
    public final CommonPrepConfig CommonPrep;
    public final String CohortPercentilesFile;

    public VisPrepConfig(final ConfigBuilder configBuilder)
    {
        CommonPrep = new CommonPrepConfig(configBuilder);
        CohortPercentilesFile = configBuilder.getValue(COHORT_PERCENTILES_FILE);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        CommonPrepConfig.registerConfig(configBuilder);
        configBuilder.addPath(COHORT_PERCENTILES_FILE, true, COHORT_PERCENTILES_FILE_DESC);
    }
}
