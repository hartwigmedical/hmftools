package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public class FindingApplicationConfig
{
    public final String OrangeJsonPath;
    public final String tumorSampleId;
    @Nullable public final String referenceSampleId;
    public final String DriverGenePath;
    public final String FindingJsonPath;
    public final String PipelineOutputDir;
    public final String PlotDir;

    private static final String ORANGE_JSON_PATH_ARG = "orange_json";
    private static final String TUMOR_SAMPLE_ID_ARG = "tumor_sample_id";
    private static final String REFERENCE_SAMPLE_ID_ARG = "reference_sample_id";
    private static final String DRIVER_GENE_PATH_ARG = "driver_gene";
    private static final String FINDING_JSON_PATH_ARG = "finding_json";
    private static final String PIPELINE_OUTPUT_DIR_ARG = "pipeline_dir";
    private static final String PLOT_DIR_ARG = "visualisation_dir";

    public FindingApplicationConfig(final ConfigBuilder configBuilder)
    {
        OrangeJsonPath =  configBuilder.getValue(ORANGE_JSON_PATH_ARG);
        tumorSampleId = configBuilder.getValue(TUMOR_SAMPLE_ID_ARG);
        referenceSampleId = configBuilder.getValue(REFERENCE_SAMPLE_ID_ARG);
        DriverGenePath =  configBuilder.getValue(DRIVER_GENE_PATH_ARG);
        FindingJsonPath = configBuilder.getValue(FINDING_JSON_PATH_ARG);
        PipelineOutputDir = configBuilder.getValue(PIPELINE_OUTPUT_DIR_ARG);
        PlotDir = configBuilder.getValue(PLOT_DIR_ARG);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(ORANGE_JSON_PATH_ARG, true, "Path to Orange JSON file");
        configBuilder.addConfigItem(TUMOR_SAMPLE_ID_ARG, true, "Tumor sample ID");
        configBuilder.addConfigItem(REFERENCE_SAMPLE_ID_ARG, false, "Reference sample ID");
        configBuilder.addPath(DRIVER_GENE_PATH_ARG, true, "Path to driver gene file");
        configBuilder.addConfigItem(FINDING_JSON_PATH_ARG, true, "Path to output finding JSON file");
        configBuilder.addPath(PIPELINE_OUTPUT_DIR_ARG, true, "Path to the pipeline output directory");
        configBuilder.addConfigItem(PLOT_DIR_ARG, true, "Path to copy plots to");

        addLoggingOptions(configBuilder);
    }
}
