package com.hartwig.hmftools.qsee.vis;

import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_DESC;
import static com.hartwig.hmftools.qsee.cohort.CohortPercentilesConfig.COHORT_PERCENTILES_FILE;
import static com.hartwig.hmftools.qsee.cohort.CohortPercentilesConfig.COHORT_PERCENTILES_FILE_DESC;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class QseePlotConfig
{
    public final String SampleFeaturesFile;

    public final String TumorId;
    public final String ReferenceId;

    public final String CohortPercentilesFile;
    public final String OutputPath;

    public static final String SAMPLE_FEATURES_FILE = "sample_features_file";
    public static final String SAMPLE_FEATURES_FILE_DESC = "Path to the sample features file";

    public static final String OUTPUT_FILE = "output_file";
    public static final String OUTPUT_FILE_DESC = "Path to write the output PDF file";

    public QseePlotConfig(final ConfigBuilder configBuilder)
    {
        SampleFeaturesFile = configBuilder.getValue(SAMPLE_FEATURES_FILE);

        TumorId = configBuilder.getValue(TUMOR);
        ReferenceId = configBuilder.getValue(REFERENCE);

        CohortPercentilesFile = configBuilder.getValue(COHORT_PERCENTILES_FILE);
        OutputPath = configBuilder.getValue(OUTPUT_FILE);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(SAMPLE_FEATURES_FILE, true, SAMPLE_FEATURES_FILE_DESC);

        configBuilder.addConfigItem(TUMOR, true, TUMOR_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);

        configBuilder.addPath(COHORT_PERCENTILES_FILE, true, COHORT_PERCENTILES_FILE_DESC);
        configBuilder.addConfigItem(OUTPUT_FILE, true, OUTPUT_FILE_DESC);

        configBuilder.addConfigItem(THREADS, false, THREADS_DESC, "1");
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
