package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public class FindingApplicationConfig
{
    public final String OrangeJsonPath;
    @Nullable public final String ClinicalTranscriptsPath;
    public final String DriverGenePath;
    public final String FindingJsonPath;

    private static final String ORANGE_JSON_PATH_ARG = "orange_json";
    private static final String CLINICAL_TRANSCRIPT_PATH_ARG = "clinical_transcript";
    private static final String DRIVER_GENE_PATH_ARG = "driver_gene";
    private static final String FINDING_JSON_PATH_ARG = "finding_json";

    public FindingApplicationConfig(final ConfigBuilder configBuilder)
    {
        OrangeJsonPath =  configBuilder.getValue(ORANGE_JSON_PATH_ARG);
        ClinicalTranscriptsPath = configBuilder.getValue(CLINICAL_TRANSCRIPT_PATH_ARG);
        DriverGenePath =  configBuilder.getValue(DRIVER_GENE_PATH_ARG);
        FindingJsonPath = configBuilder.getValue(FINDING_JSON_PATH_ARG);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(ORANGE_JSON_PATH_ARG, true, "Path to Orange JSON file");
        configBuilder.addPath(CLINICAL_TRANSCRIPT_PATH_ARG, false, "Path to clinical transcripts file");
        configBuilder.addPath(DRIVER_GENE_PATH_ARG, true, "Path to driver gene file");
        configBuilder.addConfigItem(FINDING_JSON_PATH_ARG, true, "Path to output finding JSON file");

        addLoggingOptions(configBuilder);
    }
}
