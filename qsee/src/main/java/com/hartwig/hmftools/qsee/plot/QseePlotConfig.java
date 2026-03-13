package com.hartwig.hmftools.qsee.plot;

import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.qsee.cohort.CohortPercentilesFile.COHORT_PERCENTILES_FILE_CFG;
import static com.hartwig.hmftools.qsee.cohort.CohortPercentilesFile.COHORT_PERCENTILES_FILE_CFG_DESC;

import java.util.List;

import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.qsee.common.SampleIdsLoader;

public class QseePlotConfig
{
    public final String VisDataFile;

    public final List<String> TumorIds;
    public final List<String> ReferenceIds;

    public final String CohortPercentilesFile;

    public final boolean ShowPlotWarnings;
    public final boolean MergePlots;

    public final String OutputDir;
    public final String OutputId;

    public final int Threads;

    public static final String VIS_DATA_FILE = "vis_data_file";
    public static final String VIS_DATA_FILE_DESC = "Path to the vis data file";

    public static final String SHOW_PLOT_WARNINGS = "show_plot_warnings";
    public static final String SHOW_PLOT_WARNINGS_DESC = "Show ggplot2 warnings";

    public static final String MERGE_PLOTS = "merge_plots";
    public static final String MERGE_PLOTS_DESC = "In multisample mode, merge plots into one file and delete the individual plot files";

    public QseePlotConfig(final ConfigBuilder configBuilder)
    {
        VisDataFile = configBuilder.getValue(VIS_DATA_FILE);

        SampleIdsLoader sampleIdsLoader = new SampleIdsLoader().fromConfig(configBuilder);
        TumorIds = sampleIdsLoader.tumorIds();
        ReferenceIds = sampleIdsLoader.referenceIds();

        CohortPercentilesFile = configBuilder.getValue(COHORT_PERCENTILES_FILE_CFG);

        ShowPlotWarnings = configBuilder.hasFlag(SHOW_PLOT_WARNINGS);
        MergePlots = configBuilder.hasFlag(MERGE_PLOTS);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        Threads = TaskExecutor.parseThreads(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(VIS_DATA_FILE, false, VIS_DATA_FILE_DESC);

        configBuilder.addConfigItem(TUMOR, false, TUMOR_IDS_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_IDS_DESC);
        configBuilder.addPath(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);

        configBuilder.addFlag(SHOW_PLOT_WARNINGS, SHOW_PLOT_WARNINGS_DESC);
        configBuilder.addFlag(MERGE_PLOTS, MERGE_PLOTS_DESC);

        configBuilder.addPath(COHORT_PERCENTILES_FILE_CFG, false, COHORT_PERCENTILES_FILE_CFG_DESC);
        configBuilder.addPath(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
        configBuilder.addConfigItem(OUTPUT_ID, false, OUTPUT_ID_DESC);

        configBuilder.addConfigItem(THREADS, false, THREADS_DESC, "1");
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
